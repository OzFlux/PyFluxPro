# standard Python modules
import copy
import datetime
import glob
import logging
import os
import subprocess
import tempfile
# 3rd party modules
import dateutil
import numpy
import matplotlib.pyplot as plt
# PFP modules
from scripts import constants as c
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def GapFillUsingMDS(ds, l5_info, called_by):
    """
    Purpose:
     Run the FluxNet C code to implement the MDS gap filling method.
    Usage:
    Side effects:
    Author: PRI
    Date: May 2018
    """
    l5im = l5_info[called_by]
    # get the file name
    file_path = l5_info[called_by]["info"]["file_path"]
    file_name = l5_info[called_by]["info"]["in_filename"]
    nc_full_path = os.path.join(file_path, file_name)
    nc_name = os.path.split(nc_full_path)[1]
    # get some useful metadata
    ts = int(float(ds.root["Attributes"]["time_step"]))
    site_name = ds.root["Attributes"]["site_name"]
    level = ds.root["Attributes"]["processing_level"]
    # get a temporary directory for the log, input and output filoes
    tmp_dir = tempfile.TemporaryDirectory(prefix="pfp_mds_")
    td = tmp_dir.name
    for item in ["input", "output", "log"]:
        os.makedirs(os.path.join(tmp_dir.name, item))
    # define the MDS input file location
    in_base_path = os.path.join(td, "input")
    # get a list of CSV files in the input directory
    in_path_files = glob.glob(os.path.join(in_base_path, "*.csv"))
    # and clean them out
    for in_file in in_path_files:
        if os.path.exists(in_file):
            os.remove(in_file)
    # define the MDS output file location
    out_base_path = os.path.join(td, "output", "")
    # get a list of CSV files in the output directory
    out_path_files = glob.glob(os.path.join(out_base_path, "*.csv"))
    # and clean them out
    for out_file in out_path_files:
        if os.path.exists(out_file):
            os.remove(out_file)
    # get some useful odds and ends
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    first_year = ldt["Data"][0].year
    last_year = ldt["Data"][-2].year
    # now loop over the series to be gap filled using MDS
    # open a log file for the MDS C code output
    log_file_path = os.path.join(td, "log", "mds.log")
    mdslogfile = open(log_file_path, "w")
    for fig_num, mds_label in enumerate(l5im["outputs"].keys()):
        logger.info(" Doing MDS gap filling for %s", l5im["outputs"][mds_label]["target"])
        l5im["outputs"][mds_label]["out_base_path"] = out_base_path
        l5im["outputs"][mds_label]["time_step"] = ts
        # make the output file name
        out_name = site_name+"_"+level+"_"+mds_label+"_mds.csv"
        out_file_path = os.path.join(out_base_path, out_name)
        # first, we write the yearly CSV input files
        l5im["outputs"][mds_label]["in_file_paths"] = []
        for current_year in range(first_year, last_year+1):
            in_name = nc_name.replace(".nc","_"+str(current_year)+"_MDS.csv")
            #in_name = str(current_year)+".csv"
            in_file_path = os.path.join(in_base_path, in_name)
            data, header, fmt = gfMDS_make_data_array(ds, current_year, l5im["outputs"][mds_label])
            numpy.savetxt(in_file_path, data, header=header, delimiter=",", comments="", fmt=fmt)
            l5im["outputs"][mds_label]["in_file_paths"].append(in_file_path)
        # then we construct the MDS C code command options list
        cmd = gfMDS_make_cmd_string(l5im, mds_label)
        # then we spawn a subprocess for the MDS C code
        subprocess.call(cmd, stdout=mdslogfile)
        mds_out_file = os.path.join(td, "output", "mds.csv")
        os.rename(mds_out_file, out_file_path)
        gfMDS_get_mds_output(ds, mds_label, out_file_path, l5_info, called_by)
        # mask long gaps, if requested
        gfMDS_mask_long_gaps(ds, mds_label, l5_info, called_by)
        # plot the MDS results
        target = l5im["outputs"][mds_label]["target"]
        drivers = l5im["outputs"][mds_label]["drivers"]
        title = site_name+' : Comparison of tower and MDS data for '+target
        pd = gfMDS_initplot(site_name=site_name, label=target, fig_num=fig_num,
                            title=title, nDrivers=len(drivers), show_plots=True)
        gfMDS_plot(pd, ds, mds_label, l5_info, called_by)

    # close the log file
    mdslogfile.close()
    return

def gfMDS_get_mds_output(ds, mds_label, out_file_path, l5_info, called_by):
    """
    Purpose:
     Reads the CSV file output by the MDS C code and puts the contents into
     the data structure.
    Usage:
     gfMDS_get_mds_output(ds, out_file_path, first_date, last_date, include_qc=False)
     where ds is a data structure
           out_file_path is the full path to the MDS output file
    Side effects:
     New series are created in the data structure to hold the MDS data.
    Author: PRI
    Date: May 2018
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # get the MDS flag value from the processing level
    processing_level = ds.root["Attributes"]["processing_level"]
    level_number = pfp_utils.strip_non_numeric(processing_level)
    # MDS flag will be 470 at L4, 570 and L5
    mds_flag_value = int(level_number)*100 + 70
    # get the name for the description variable attribute
    descr_level = "description_" + str(ds.root["Attributes"]["processing_level"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    data_mds = numpy.genfromtxt(out_file_path, delimiter=",", names=True, autostrip=True, dtype=None)
    dt_mds = numpy.array([dateutil.parser.parse(str(dt)) for dt in data_mds["TIMESTAMP"]])
    idxa, idxb = pfp_utils.FindMatchingIndices(ldt["Data"], dt_mds)
    # get a list of the names in the data array
    mds_output_names = list(data_mds.dtype.names)
    # strip out the timestamp and the original data
    for item in ["TIMESTAMP", l5_info[called_by]["outputs"][mds_label]["target_mds"]]:
        if item in mds_output_names:
            mds_output_names.remove(item)
    # and now loop over the MDS output series
    for mds_output_name in mds_output_names:
        if mds_output_name == "FILLED":
            # get the gap filled target and write it to the data structure
            var_in = pfp_utils.GetVariable(ds, l5_info[called_by]["outputs"][mds_label]["target"])
            data_out = numpy.ma.array(var_in["Data"], copy=True)
            data_out[idxa] = data_mds[mds_output_name][idxb]
            data_out = numpy.ma.masked_values(data_out, c.missing_value)
            # Ugly hack for datasets that start on YYYY-01-01 00:00, the MDS C code discards
            # these times and they are missing from the MDS output.
            # Many thanks to Professor Shih-Chieh Chang of National Dong Hwa University, Taiwan
            # for reporting this bug.  Problems often behave like fractals.
            if numpy.ma.is_masked(data_out[0]):
                data_out[0] = data_out[1]
            idx = numpy.where(numpy.ma.getmaskarray(var_in["Data"]) == True)[0]
            flag = numpy.array(var_in["Flag"], copy=True)
            flag[idx] = numpy.int32(mds_flag_value)
            attr = copy.deepcopy(var_in["Attr"])
            pfp_utils.append_to_attribute(attr, {descr_level: "Gap filled using MDS"})
            var_out = {"Label": mds_label, "Data": data_out, "Flag": flag, "Attr": attr}
            pfp_utils.CreateVariable(ds, var_out)
        elif mds_output_name == "TIMEWINDOW":
            # make the series name for the data structure
            mds_qc_label = "MDS"+"_"+l5_info[called_by]["outputs"][mds_label]["target"]+"_"+mds_output_name
            data_out = numpy.zeros(nrecs, dtype=float)
            data_out[idxa] = data_mds[mds_output_name][idxb]
            flag = numpy.zeros(nrecs, dtype=int)
            attr = {"long_name":"TIMEWINDOW from MDS gap filling for "+l5_info[called_by]["outputs"][mds_label]["target"]}
            var_out = {"Label": mds_qc_label, "Data": data_out, "Flag": flag, "Attr": attr}
            pfp_utils.CreateVariable(ds, var_out)
        else:
            # make the series name for the data structure
            mds_qc_label = "MDS"+"_"+l5_info[called_by]["outputs"][mds_label]["target"]+"_"+mds_output_name
            data_out = numpy.zeros(nrecs, dtype=float)
            data_out[idxa] = data_mds[mds_output_name][idxb]
            flag = numpy.zeros(nrecs, dtype=int)
            attr = {"long_name":"QC field from MDS gap filling for "+l5_info[called_by]["outputs"][mds_label]["target"]}
            var_out = {"Label": mds_qc_label, "Data": data_out, "Flag": flag, "Attr": attr}
            pfp_utils.CreateVariable(ds, var_out)
    return

def gfMDS_initplot(**kwargs):
    # set the margins, heights, widths etc
    pd = {"margin_bottom":0.075, "margin_top":0.075, "margin_left":0.05, "margin_right":0.05,
          "xy_height":0.20, "xy_width":0.20, "xyts_space":0.05, "ts_width":0.9}
    # set the keyword arguments
    for key, value in kwargs.items():
        pd[key] = value
    # calculate bottom of the first time series and the height of the time series plots
    pd["ts_bottom"] = pd["margin_bottom"]+pd["xy_height"]+pd["xyts_space"]
    pd["ts_height"] = (1.0 - pd["margin_top"] - pd["ts_bottom"])/float(pd["nDrivers"]+1)
    return pd

def gfMDS_make_cmd_string(l5im, mds_label):
    """
    Purpose:
     Construct the command line options string passed to the MDS C code.
    Usage:
    Author: PRI
    Date: May 2018
    """
    info = l5im["outputs"][mds_label]
    # create the input files string
    in_files = info["in_file_paths"][0]
    suffix = l5im["info"]["executable_suffix"]
    for in_file in info["in_file_paths"][1:]:
        in_files = in_files + "+" + in_file
    # get the output base path
    out_base_path = info["out_base_path"]
    # get the base path of script or Pyinstaller application
    base_path = pfp_utils.get_base_path()
    gf_mds_exe = os.path.join(base_path, "mds", "bin", "gf_mds"+suffix)
    # start creating the command list of MDS options
    cmd = [gf_mds_exe, '-input='+in_files, '-output='+out_base_path,
           '-date=TIMESTAMP', '-rows_min=0']
    # create the target label string
    tofill = info["target_mds"]
    cmd.append('-tofill='+tofill)
    # create the driver label and tolerance strings
    sw_in = info["drivers_mds"][0]
    cmd.append('-sw_in='+sw_in)
    sw_int = str(info["tolerances"][0][0]) + ","
    sw_int = sw_int + str(info["tolerances"][0][1])
    cmd.append('-sw_int='+sw_int)
    if len(info["drivers_mds"][1]) > 0:
        ta = info["drivers_mds"][1]
        cmd.append('-ta='+ta)
        tat = str(info["tolerances"][1])
        cmd.append('-tat='+tat)
    if len(info["drivers_mds"][2]):
        vpd = info["drivers_mds"][2]
        cmd.append('-vpd='+vpd)
        vpdt = str(info["tolerances"][2])
        cmd.append('-vpdt='+vpdt)
    if info["time_step"] == 60:
        cmd.append('-hourly')
    return cmd

def gfMDS_make_data_array(ds, current_year, info):
    """
    Purpose:
     Create a data array for the MDS gap filling routine.  The array constructed
     here will be written to a CSV file that is read by the MDS C code.
    Usage:
    Side Effects:
     The constructed data arrays are full years.  That is they run from YYYY-01-01 00:30
     to YYYY+1-01-01 00:00.  Missing data is represented as -9999.
    Author: PRI
    Date: May 2018
    """
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ts = int(float(ds.root["Attributes"]["time_step"]))
    start = datetime.datetime(current_year, 1, 1, 0, 0, 0) + datetime.timedelta(minutes=ts)
    end = datetime.datetime(current_year+1, 1, 1, 0, 0, 0)
    cdt = numpy.array([dt for dt in pfp_utils.perdelta(start, end, datetime.timedelta(minutes=ts))])
    mt = numpy.ones(len(cdt))*float(-9999)
    # need entry for the timestamp and the target ...
    array_list = [cdt, mt]
    # ... and entries for the drivers
    for driver in info["drivers"]:
        array_list.append(mt)
    # now we can create the data array
    data = numpy.stack(array_list, axis=-1)
    si = pfp_utils.GetDateIndex(ldt["Data"], start, default=0)
    ei = pfp_utils.GetDateIndex(ldt["Data"], end, default=nrecs)
    dt = pfp_utils.GetVariable(ds, "DateTime", start=si, end=ei)
    idx1, _ = pfp_utils.FindMatchingIndices(cdt, dt["Data"])
    pfp_label_list = [info["target"]]+info["drivers"]
    mds_label_list = [info["target_mds"]]+info["drivers_mds"]
    header = "TIMESTAMP"
    fmt = "%12i"
    for n, label in enumerate(pfp_label_list):
        var = pfp_utils.GetVariable(ds, label, start=si, end=ei)
        data[idx1,n+1] = var["Data"]
        header = header+","+mds_label_list[n]
        fmt = fmt+","+"%f"
    # convert datetime to ISO dates
    data[:,0] = numpy.array([int(xdt.strftime("%Y%m%d%H%M")) for xdt in cdt])
    return data, header, fmt

def gfMDS_mask_long_gaps(ds, mds_label, l5_info, called_by):
    """
    Purpose:
     Mask gaps that are longer than a specified maximum length.
    Usage:
    Side effects:
    Author: PRI
    Date: June 2019
    """
    if not l5_info[called_by]["outputs"][mds_label]["mask long gaps"]:
        # mask long gaps disabled in control file
        return
    if "MaxShortGapRecords" not in l5_info[called_by]["info"]:
        return
    max_short_gap_days = l5_info[called_by]["info"]["MaxShortGapDays"]
    msg = "  Masking gaps longer than " + str(max_short_gap_days) + " days"
    logger.info(msg)
    label = l5_info[called_by]["outputs"][mds_label]["target"]
    target = pfp_utils.GetVariable(ds, label)
    variable = pfp_utils.GetVariable(ds, mds_label)
    mask = numpy.ma.getmaskarray(target["Data"])
    # start and stop indices of contiguous blocks
    max_short_gap_records = l5_info[called_by]["info"]["MaxShortGapRecords"]
    gap_start_end = pfp_utils.contiguous_regions(mask)
    for start, stop in gap_start_end:
        gap_length = stop - start
        if gap_length > max_short_gap_records:
            variable["Data"][start: stop] = target["Data"][start: stop]
            variable["Flag"][start: stop] = target["Flag"][start: stop]
    # put data_int back into the data structure
    pfp_utils.CreateVariable(ds, variable)
    return

def gfMDS_plot(pd, ds, mds_label, l5_info, called_by):
    """
    Purpose:
     Plot the drivers, target and gap filled variable.
    Usage:
    Author: PRI
    Date: Back in the day
    """
    # get the MDS flag value from the processing level
    processing_level = ds.root["Attributes"]["processing_level"]
    level_number = pfp_utils.strip_non_numeric(processing_level)
    # MDS flag will be 470 at L4, 570 and L5
    mds_flag_value = int(level_number)*100 + 70
    # get the timestep
    ts = int(float(ds.root["Attributes"]["time_step"]))
    # get a local copy of the datetime series
    xdt = ds.root["Variables"]["DateTime"]["Data"]
    Hdh = numpy.array([d.hour+(d.minute+d.second/float(60))/float(60) for d in xdt])
    # get the observed and modelled values
    drivers = l5_info[called_by]["outputs"][mds_label]["drivers"]
    target = l5_info[called_by]["outputs"][mds_label]["target"]
    obs = pfp_utils.GetVariable(ds, target)
    mds = pfp_utils.GetVariable(ds, mds_label)
    if pd["show_plots"]:
        plt.ion()
    else:
        plt.ioff()
    if plt.fignum_exists(1):
        fig = plt.figure(1)
        plt.clf()
    else:
        fig = plt.figure(1, figsize=(13, 8))
    fig.canvas.manager.set_window_title(target)
    plt.figtext(0.5, 0.95, pd["title"], ha='center', size=16)

    # diurnal plot
    # XY plot of the diurnal variation
    rect1 = [0.10, pd["margin_bottom"], pd["xy_width"], pd["xy_height"]]
    ax1 = plt.axes(rect1)
    # get the diurnal stats of the observations
    mask = numpy.ma.mask_or(obs["Data"].mask, mds["Data"].mask)
    obs_mor = numpy.ma.array(obs["Data"], mask=mask, copy=True)
    _, Hr1, Av1, _, _, _ = gf_getdiurnalstats(Hdh, obs_mor, ts)
    ax1.plot(Hr1, Av1, 'b-', label="Obs")
    # get the diurnal stats of all SOLO predictions
    _, Hr2, Av2, _, _, _ = gf_getdiurnalstats(Hdh, mds["Data"], ts)
    ax1.plot(Hr2, Av2, 'r-', label="MDS")
    plt.xlim(0, 24)
    plt.xticks([0, 6, 12, 18, 24])
    ax1.set_ylabel(target)
    ax1.set_xlabel('Hour')
    ax1.legend(loc='upper right', frameon=False, prop={'size':8})

    # histogram of window size
    time_window = pfp_utils.GetVariable(ds, "MDS_"+target+"_TIMEWINDOW")
    idx = numpy.where(mds["Flag"] == mds_flag_value)[0]
    if len(idx) != 0:
        tw_hist_data = time_window["Data"][idx]
        rect2 = [0.40,pd["margin_bottom"],pd["xy_width"],pd["xy_height"]]
        ax2 = plt.axes(rect2)
        ax2.hist(tw_hist_data)
        ax2.set_ylabel("Occurrence")
        ax2.set_xlabel("MDS window length")

    # write statistics to the plot
    numpoints = numpy.ma.count(obs["Data"])
    numfilled = numpy.ma.count(mds["Data"])-numpy.ma.count(obs["Data"])
    plt.figtext(0.65,0.225,'No. points')
    plt.figtext(0.75,0.225,str(numpoints))
    plt.figtext(0.65,0.200,'No. filled')
    plt.figtext(0.75,0.200,str(numfilled))
    avg_obs = numpy.ma.mean(obs["Data"])
    avg_mds = numpy.ma.mean(mds["Data"])
    plt.figtext(0.65,0.175,'Avg (obs)')
    plt.figtext(0.75,0.175,'%.4g'%(avg_obs))
    plt.figtext(0.65,0.150,'Avg (MDS)')
    plt.figtext(0.75,0.150,'%.4g'%(avg_mds))
    var_obs = numpy.ma.var(obs["Data"])
    var_mds = numpy.ma.var(mds["Data"])
    plt.figtext(0.65,0.125,'Var (obs)')
    plt.figtext(0.75,0.125,'%.4g'%(var_obs))
    plt.figtext(0.65,0.100,'Var (MDS)')
    plt.figtext(0.75,0.100,'%.4g'%(var_mds))

    # time series of drivers and target
    ts_axes = []
    rect = [pd["margin_left"],pd["ts_bottom"],pd["ts_width"],pd["ts_height"]]
    ts_axes.append(plt.axes(rect))
    ts_axes[0].plot(obs["DateTime"], obs["Data"], 'b.',
                    mds["DateTime"], mds["Data"], 'r-')
    ts_axes[0].set_xlim(obs["DateTime"][0], obs["DateTime"][-1])
    TextStr = target+'_obs ('+obs['Attr']['units']+')'
    ts_axes[0].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[0].transAxes)
    TextStr = target+'('+mds['Attr']['units']+')'
    ts_axes[0].text(0.85,0.85,TextStr,color='r',horizontalalignment='right',transform=ts_axes[0].transAxes)
    for i, driver in enumerate(drivers):
        this_bottom = pd["ts_bottom"] + (i+1)*pd["ts_height"]
        rect = [pd["margin_left"], this_bottom, pd["ts_width"], pd["ts_height"]]
        ts_axes.append(plt.axes(rect, sharex=ts_axes[0]))
        drv = pfp_utils.GetVariable(ds, driver)
        drv_notgf = numpy.ma.masked_where(drv["Flag"] != 0, drv["Data"])
        drv_gf = numpy.ma.masked_where(drv["Flag"] == 0, drv["Data"])
        ts_axes[i+1].plot(drv["DateTime"], drv_notgf, 'b-')
        ts_axes[i+1].plot(drv["DateTime"], drv_gf, 'r-', linewidth=2)
        plt.setp(ts_axes[i+1].get_xticklabels(), visible=False)
        TextStr = driver+'('+drv['Attr']['units']+')'
        ts_axes[i+1].text(0.05,0.85,TextStr,color='b',horizontalalignment='left',transform=ts_axes[i+1].transAxes)

    # save a hard copy
    sdt = obs["DateTime"][0].strftime("%Y%m%d")
    edt = obs["DateTime"][-1].strftime("%Y%m%d")
    figure_name = pd["site_name"].replace(" ","") + "_MDS_" + pd["label"]
    figure_name = figure_name + "_" + sdt + "_" + edt + ".png"
    figure_path = os.path.join(l5_info[called_by]["info"]["plot_path"], figure_name)
    fig.savefig(figure_path, format='png')
    if pd["show_plots"]:
        fig.canvas.flush_events()
    else:
        plt.close(fig)
        plt.ion()
    return

def gf_getdiurnalstats(DecHour, Data, ts):
    nInts = 24*int((60/ts)+0.5)
    Num = numpy.ma.zeros(nInts, dtype=int)
    Hr = numpy.ma.zeros(nInts, dtype=float)
    for i in range(nInts):
        Hr[i] = float(i)*ts/60.
    Av = numpy.ma.masked_all(nInts)
    Sd = numpy.ma.masked_all(nInts)
    Mx = numpy.ma.masked_all(nInts)
    Mn = numpy.ma.masked_all(nInts)
    if numpy.size(Data)!=0:
        for i in range(nInts):
            li = numpy.ma.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
            Num[i] = numpy.size(li)
            if Num[i]!=0:
                Av[i] = numpy.ma.mean(Data[li])
                Sd[i] = numpy.ma.std(Data[li])
                Mx[i] = numpy.ma.max(Data[li])
                Mn[i] = numpy.ma.min(Data[li])
    return Num, Hr, Av, Sd, Mx, Mn
