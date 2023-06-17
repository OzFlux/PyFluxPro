# Python modules
import copy
import datetime
import logging
import os
# 3rd party modules
import matplotlib.pyplot as plt
import numpy
import scipy
import xlsxwriter
# PFP modules
from scripts import pfp_classes
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def check_final_e0(results, year):
    """ Based on ONEFlux code for estimating final rb values."""
    rwqc = results["years"][year]["winqc"]
    not_nan = numpy.where(~numpy.isnan(rwqc["E0"]["value"]) &
                          ~numpy.isnan(rwqc["E0"]["stderr"]))[0]
    e0_nn = {"value": rwqc["E0"]["value"][not_nan],
             "stderr": rwqc["E0"]["stderr"][not_nan]}
    idx = numpy.where((e0_nn["stderr"] < 100) &
                      (e0_nn["stderr"]/e0_nn["value"] < 0.5) &
                      (e0_nn["value"] > 50) &
                      (e0_nn["stderr"] < 450))[0]
    if len(idx) == 0:
        msg = " No E0 from windows passed quality control, skipping year " + str(year)
        logger.warning(msg)
        ok = False
    else:
        ok = True
    return ok
def create_results(years):
    """ Create a dictionary to hold the results. """
    results = {"years": {}, "final": {}}
    mt = numpy.array([])
    for year in years:
        results["years"][year] = {}
        results["years"][year]["all"] = {"E0": {"value": mt, "stderr": mt},
                                         "rb": {"value": mt, "stderr": mt}}
        results["years"][year]["allqc"] = {"E0": {"value": mt, "stderr": mt},
                                           "rb": {"value": mt, "stderr": mt}}
        results["years"][year]["win"] = {"start": mt, "mid": mt, "end": mt,
                                         "num": mt, "Trange": mt,
                                         "E0": {"value": mt, "stderr": mt},
                                         "rb": {"value": mt, "stderr": mt}}
        results["years"][year]["winqc"] = {"start": mt, "mid": mt, "end": mt,
                                           "num": mt, "Trange": mt,
                                           "E0": {"value": mt, "stderr": mt},
                                           "rb": {"value": mt, "stderr": mt}}
        results["years"][year]["winqc_avg"] = {"E0": {"value": mt, "stderr": mt}}
        results["years"][year]["final"] = {"start": mt, "mid": mt, "end": mt,
                                           "num": mt, "Trange": mt,
                                           "E0": {"value": mt, "stderr": mt},
                                           "rb": {"value": mt, "stderr": mt}}
    results["final"] = {"start": mt, "mid": mt, "end": mt,
                        "num": mt, "Trange": mt,
                        "E0": {"value": mt, "stderr": mt},
                        "rb": {"value": mt, "stderr": mt}}
    results["win"] = {"start": mt, "mid": mt, "end": mt,
                      "num": mt, "Trange": mt,
                      "E0": {"value": mt, "stderr": mt},
                      "rb": {"value": mt, "stderr": mt}}
    return results
def er_lloyd_taylor(T, E0, rb):
    return rb*numpy.exp(E0*(1/(288.13-227.13)-1/((T+273.15)-227.13)))
def estimate_e0_full_year(ds_year, results, l6_info):
    called_by = l6_info["ERUsingLloydTaylor"]["info"]["called_by"]
    ts = int(ds_year.root["Attributes"]["time_step"])
    ER = pfp_utils.GetVariable(ds_year, "ER")
    Ta = pfp_utils.GetVariable(ds_year, "Ta")
    ldt = ds_year.root["Variables"]["DateTime"]["Data"] - datetime.timedelta(minutes=ts)
    year = numpy.unique([dt.year for dt in ldt])[0]
    e0p = int(l6_info[called_by]["outputs"]["ER_LT_all"]["e0_prior"])
    rbp = int(l6_info[called_by]["outputs"]["ER_LT_all"]["rb_prior"])
    mask = numpy.ma.getmaskarray(ER["Data"])
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Ta["Data"]), copy=True, shrink=False)
    y = numpy.ma.compressed(numpy.ma.masked_where(mask, ER["Data"]))
    x = numpy.ma.compressed(numpy.ma.masked_where(mask, Ta["Data"]))
    rls = least_squares(residuals_lt, numpy.array([e0p, rbp]), (x, y), method='lm')
    res = results["years"][year]["all"]
    res["E0"]["value"] = numpy.append(res["E0"]["value"], rls["rls"].x[0])
    res["E0"]["stderr"] = numpy.append(res["E0"]["stderr"], rls["perr"][0])
    res["rb"]["value"] = numpy.append(res["rb"]["value"], rls["rls"].x[1])
    res["rb"]["stderr"] = numpy.append(res["rb"]["stderr"], rls["perr"][1])
    return
def estimate_e0_windows(ds_year, results, l6_info, fit_type="least_squares"):
    called_by = l6_info["ERUsingLloydTaylor"]["info"]["called_by"]
    nrecs = int(ds_year.root["Attributes"]["nc_nrecs"])
    ts = int(ds_year.root["Attributes"]["time_step"])
    nperday = int(24*60/ts)
    ldt = ds_year.root["Variables"]["DateTime"]["Data"] - datetime.timedelta(minutes=ts)
    year = numpy.unique([dt.year for dt in ldt])[0]
    ws = int(l6_info[called_by]["outputs"]["ER_LT_all"]["e0_window_size_days"])
    ss = int(l6_info[called_by]["outputs"]["ER_LT_all"]["e0_step_size_days"])
    mts = int(l6_info[called_by]["outputs"]["ER_LT_all"]["minimum_temperature_spread"])
    mwp = int(l6_info[called_by]["outputs"]["ER_LT_all"]["minimum_window_points"])
    e0p = int(l6_info[called_by]["outputs"]["ER_LT_all"]["e0_prior"])
    rbp = int(l6_info[called_by]["outputs"]["ER_LT_all"]["rb_prior"])
    ER = pfp_utils.GetVariable(ds_year, "ER")
    Ta = pfp_utils.GetVariable(ds_year, "Ta")
    Sws = pfp_utils.GetVariable(ds_year, "Sws")
    mask = numpy.ma.getmaskarray(ER["Data"])
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Ta["Data"]))
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Sws["Data"]))
    ER["Data"] = numpy.ma.masked_where(mask, ER["Data"])
    Ta["Data"] = numpy.ma.masked_where(mask, Ta["Data"])
    Sws["Data"] = numpy.ma.masked_where(mask, Sws["Data"])
    si = 0
    ei = si + ws*nperday
    n = 0
    while ei < nrecs:
        y = numpy.ma.compressed(ER["Data"][si:ei+1])
        x = numpy.ma.compressed(Ta["Data"][si:ei+1])
        z = numpy.ma.compressed(Sws["Data"][si:ei+1])
        if (len(y) < mwp):
            msg = "  Insufficient number of points " + str(len(y))
            logger.debug(msg)
            continue
        if max(x) - min(x) < mts:
            msg = "  Insufficient temperature spread " + str(max(y)-min(y))
            logger.debug(msg)
            continue
        rls = least_squares(residuals_lt, numpy.array([e0p, rbp]), (x, y), method='lm')
        rw = results["years"][year]["win"]
        rw["start"] = numpy.append(rw["start"], ER["DateTime"][si])
        rw["end"] = numpy.append(rw["end"], ER["DateTime"][ei])
        mid = rw["start"][-1] + (rw["end"][-1]-rw["start"][-1])/2
        rw["mid"] = numpy.append(rw["mid"], mid)
        rw["num"] = numpy.append(rw["num"], n)
        rw["Trange"] = numpy.append(rw["Trange"], max(x) - min(x))
        rw["E0"]["value"] = numpy.append(rw["E0"]["value"], rls["rls"].x[0])
        rw["E0"]["stderr"] = numpy.append(rw["E0"]["stderr"], rls["perr"][0])
        rw["rb"]["value"] = numpy.append(rw["rb"]["value"], rls["rls"].x[1])
        rw["rb"]["stderr"] = numpy.append(rw["rb"]["stderr"], rls["perr"][1])
        # plot to screen and hard copy
        if l6_info["Options"]["plot_raw_data"]:
            plot_er_vs_ta(y, x, z, rw, l6_info)
        si = si + ss*nperday
        ei = si + ws*nperday
        n = n + 1
    return
def estimate_er_lt(ds, output):
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ER_LT = pfp_utils.CreateEmptyVariable(output, nrecs)
    Ta = pfp_utils.GetVariable(ds, "Ta")
    E0 = pfp_utils.GetVariable(ds, "E0_LT")
    Rb = pfp_utils.GetVariable(ds, "Rb_LT")
    ER_LT["Data"] = er_lloyd_taylor(Ta["Data"], E0["Data"], Rb["Data"])
    ER_LT["Attr"]["long_name"] = "Ecosystem respiration"
    ER_LT["Attr"]["description_L6"] = "Ecosystem respiration modelled by Lloyd-Taylor"
    ER_LT["Attr"]["units"] = "umol/m^2/s"
    ER_LT["Attr"]["statistic_type"] = "average"
    pfp_utils.CreateVariable(ds, ER_LT)
    return
def get_final_e0(results, year):
    """ Based on ONEFlux code for selecting final E0 value."""
    raqc = results["years"][year]["allqc"]
    rwqc = results["years"][year]["winqc"]
    rwqca = results["years"][year]["winqc_avg"]
    idx = ~numpy.isnan(rwqc["E0"]["value"])
    e0_not_nan = rwqc["E0"]["value"][idx]
    if len(e0_not_nan) > 1:
        e0_value = rwqc["E0"]["value"][idx]
        e0_stderr = rwqc["E0"]["stderr"][idx]
        e0_results = numpy.column_stack([e0_value, e0_stderr])
        e0_sorted = e0_results[numpy.argsort(e0_results[:, -1])]
        n = min([3, e0_sorted.shape[0]])
        rwqca["E0"]["value"] = numpy.mean(e0_sorted[:n, 0])
        rwqca["E0"]["stderr"] = numpy.mean(e0_sorted[:n, 1])
    else:
        rwqca["E0"]["value"] = raqc["E0"]["value"]
        rwqca["E0"]["stderr"] = raqc["E0"]["stderr"]
    return
def get_final_rb(ds_year, results, l6_info):
    called_by = l6_info["ERUsingLloydTaylor"]["info"]["called_by"]
    nrecs = int(ds_year.root["Attributes"]["nc_nrecs"])
    ts = int(ds_year.root["Attributes"]["time_step"])
    nperday = int(24*60/ts)
    ldt = ds_year.root["Variables"]["DateTime"]["Data"] - datetime.timedelta(minutes=ts)
    year = numpy.unique([dt.year for dt in ldt])[0]
    e0 = results["years"][year]["winqc_avg"]["E0"]["value"]
    ws = int(l6_info[called_by]["outputs"]["ER_LT_all"]["rb_window_size_days"])
    ss = int(l6_info[called_by]["outputs"]["ER_LT_all"]["rb_step_size_days"])
    mts = int(l6_info[called_by]["outputs"]["ER_LT_all"]["minimum_temperature_spread"])
    rbp = int(l6_info[called_by]["outputs"]["ER_LT_all"]["rb_prior"])
    ER = pfp_utils.GetVariable(ds_year, "ER")
    Ta = pfp_utils.GetVariable(ds_year, "Ta")
    mask = numpy.ma.getmaskarray(ER["Data"])
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Ta["Data"]))
    ER["Data"] = numpy.ma.masked_where(mask, ER["Data"])
    Ta["Data"] = numpy.ma.masked_where(mask, Ta["Data"])
    si = 0
    ei = si + ws*nperday
    n = 0
    while ei < nrecs:
        y = numpy.ma.compressed(ER["Data"][si:ei+1])
        x = numpy.ma.compressed(Ta["Data"][si:ei+1])
        if ((len(y) >= 6) and (max(x) - min(x) >= mts)):
            rls = least_squares(residuals_lt_rb, numpy.array([rbp]), (x, y, e0), method='lm')
            rf = results["years"][year]["final"]
            rf["start"] = numpy.append(rf["start"], ER["DateTime"][si])
            rf["end"] = numpy.append(rf["end"], ER["DateTime"][ei])
            mid = rf["start"][-1] + (rf["end"][-1]-rf["start"][-1])/2
            rf["mid"] = numpy.append(rf["mid"], mid)
            rf["num"] = numpy.append(rf["num"], len(y))
            rf["Trange"] = numpy.append(rf["Trange"], max(x) - min(x))
            rf["E0"]["value"] = numpy.append(rf["E0"]["value"], e0)
            rf["E0"]["stderr"] = numpy.append(rf["E0"]["stderr"], float(0))
            rf["rb"]["value"] = numpy.append(rf["rb"]["value"], rls["rls"].x[0])
            rf["rb"]["stderr"] = numpy.append(rf["rb"]["stderr"], rls["perr"][0])
        elif (len(y) < 6):
            msg = "  Insufficient number of points " + str(len(y))
            logger.debug(msg)
        elif (max(x) - min(x) < mts):
            msg = "  Insufficient temperature spread " + str(max(y)-min(y))
            logger.debug(msg)
        # plot to screen and hard copy
        # plot_er_vs_ta(y, x, z, l6_info)
        si = si + ss*nperday
        ei = si + ws*nperday
        n = n + 1
    return
def interpolate_1d(ts_out, ts_in, data_in, kind='linear', bounds_error=False, fill_value=numpy.nan):
    """ Duplicates the ONEFlux nighttime.ipolmiss() function. """
    not_nan = numpy.where(~numpy.isnan(data_in))[0]
    ifunc = scipy.interpolate.interp1d(ts_in[not_nan], data_in[not_nan],
                                       kind=kind, bounds_error=bounds_error, fill_value=fill_value)
    data_out = ifunc(ts_out)
    # extend first value to beginning to time series and last value to end
    not_nan = numpy.where(~numpy.isnan(data_out))[0]
    data_out[0:not_nan[0]] = data_out[not_nan[0]]
    data_out[not_nan[-1]+1:-1] = data_out[not_nan[-1]]
    return data_out
def interpolate_parameters_lt(ds, results):
    """ Interpolate the LT parameters from the window time step onto the tower time step. """
    # get a list of years
    years = list(results["years"].keys())
    # loop over years and concatenate results into single array for each parameter
    rf = results["final"]
    for year in years:
        ryf = results["years"][year]["final"]
        rf["start"] = numpy.append(rf["start"], ryf["start"])
        rf["mid"] = numpy.append(rf["mid"], ryf["mid"])
        rf["end"] = numpy.append(rf["end"], ryf["end"])
        rf["E0"]["value"] = numpy.append(rf["E0"]["value"], ryf["E0"]["value"])
        rf["E0"]["stderr"] = numpy.append(rf["E0"]["stderr"], ryf["E0"]["stderr"])
        rf["rb"]["value"] = numpy.append(rf["rb"]["value"], ryf["rb"]["value"])
        rf["rb"]["stderr"] = numpy.append(rf["rb"]["stderr"], ryf["rb"]["stderr"])
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # create the variables for the interpolated parameters
    attr = {"long_name": "Basal respiration", "units": "umol/m^2/s"}
    Rb = pfp_utils.CreateEmptyVariable("Rb_LT", nrecs, attr=attr)
    attr = {"long_name": "Basal respiration", "statistic_type": "standard_deviation", "units": "umol/m^2/s"}
    Rb_se = pfp_utils.CreateEmptyVariable("Rb_LT_se", nrecs, attr=attr)
    attr = {"long_name": "Activation energy", "units": "K^-1"}
    E0 = pfp_utils.CreateEmptyVariable("E0_LT", nrecs, attr=attr)
    attr = {"long_name": "Activation energy", "statistic_type": "standard_deviation", "units": "K^-1"}
    E0_se = pfp_utils.CreateEmptyVariable("E0_LT_se", nrecs, attr=attr)
    # get the tower time step and the mid-window time step
    ts_ds = numpy.array([d.timestamp() for d in ds.root["Variables"]["DateTime"]["Data"]])
    ts_ltr = numpy.array([d.timestamp() for d in rf["mid"]])
    # interpolate Rb value onto the tower time step
    Rb["Data"] = interpolate_1d(ts_ds, ts_ltr, rf["rb"]["value"])
    pfp_utils.CreateVariable(ds, Rb)
    # interpolate Rb standard error onto the tower time step
    Rb_se["Data"] = interpolate_1d(ts_ds, ts_ltr, rf["rb"]["stderr"])
    pfp_utils.CreateVariable(ds, Rb_se)
    # interpolate E0 value onto the tower time step
    E0["Data"] = interpolate_1d(ts_ds, ts_ltr, rf["E0"]["value"])
    pfp_utils.CreateVariable(ds, E0)
    # interpolate E0 standard error onto the tower time step
    E0_se["Data"] = interpolate_1d(ts_ds, ts_ltr, rf["E0"]["stderr"])
    pfp_utils.CreateVariable(ds, E0_se)
    return
def least_squares(func, p0, args, method='trf', loss='linear'):
    rls = scipy.optimize.least_squares(func, p0, args=args, method=method, loss=loss)
    # see https://stackoverflow.com/questions/42388139
    U, s, Vh = numpy.linalg.svd(rls.jac, full_matrices=False)
    tol = numpy.finfo(float).eps*s[0]*max(rls.jac.shape)
    w = s > tol
    # robust covariance matrix
    cov = (Vh[w].T/s[w]**2) @ Vh[w]
    # 1sigma uncertainty on fitted parameters
    perr = numpy.sqrt(numpy.diag(cov))
    return {"rls": rls, "perr": perr}
def plot_er_vs_ta(er, ta, sws, rw, l6_info):
    fignum = l6_info["ERUsingLloydTaylor"]["info"]["fignum"]
    title = str(rw["start"][-1]) + " to " + str(rw["end"][-1])
    plot_path = l6_info["ERUsingLloydTaylor"]["info"]["plot_path"]
    if not os.path.isdir(plot_path):
        os.mkdir(plot_path)
    file_name = os.path.join(plot_path, "estimate_e0_" + title.replace(" ", "_") + ".png")
    if l6_info["Options"]["call_mode"] == "interactive":
        plt.ion()
    else:
        current_backend = plt.get_backend()
        plt.switch_backend("agg")
        plt.ioff()
    if plt.fignum_exists(fignum):
        fig = plt.figure(fignum)
        plt.clf()
        axs = fig.add_subplot(1, 1, 1)
    else:
        fig, axs = plt.subplots(num=fignum, figsize=(6, 6))
    sc = axs.scatter(ta, er, c=sws, s=10)
    ta_fit = numpy.linspace(min(ta), max(ta), 10)
    er_fit = er_lloyd_taylor(ta_fit, rw["E0"]["value"][-1], rw["rb"]["value"][-1])
    axs.plot(ta_fit, er_fit, 'r--')
    axs.set_title(title)
    axs.set_xlabel("Temperature (degC)")
    axs.set_ylabel("ER (umol/m^2/s)")
    clb = plt.colorbar(sc)
    clb.ax.set_title("Sws")
    fig.savefig(file_name, format="png")
    if l6_info["Options"]["call_mode"] == "interactive":
        plt.draw()
        pfp_utils.mypause(0.5)
        plt.ioff()
    else:
        plt.close()
        plt.switch_backend(current_backend)
        plt.ion()
    return
def quality_control_e0(results, year):
    """ Based on ONEFlux quality control of E0."""
    # qc values based on whole year of data
    ra = results["years"][year]["all"]
    raqc = results["years"][year]["allqc"] = copy.deepcopy(ra)
    # ONEFlux limits 0 < E0 < 450 for E0 values based on whole year data
    idx = numpy.where(ra["E0"]["value"] < 0)[0]
    raqc["E0"]["value"][idx] = float(0)
    raqc["E0"]["stderr"][idx] = float(0)
    idx = numpy.where(ra["E0"]["value"] > 450)[0]
    raqc["E0"]["value"][idx] = float(450)
    raqc["E0"]["stderr"][idx] = float(0)
    # qc values based on windowed data
    rw = results["years"][year]["win"]
    rwqc = results["years"][year]["winqc"] = copy.deepcopy(rw)
    # ONEFlux limits 30 < E0 < 450 for E0 values based on windowed data
    idx = numpy.where((rw["E0"]["value"] < 30) | (rw["E0"]["value"] > 350))[0]
    rwqc["E0"]["value"][idx] = numpy.nan
    rwqc["E0"]["stderr"][idx] = numpy.nan
    rwqc["rb"]["value"][idx] = numpy.nan
    rwqc["rb"]["stderr"][idx] = numpy.nan
    return
def residuals_lt(params, T, ER):
    E0 = params[0]
    rb = params[1]
    residuals = er_lloyd_taylor(T, E0, rb) - ER
    return residuals
def residuals_lt_rb(params, T, ER, E0):
    rb = params[0]
    residuals = er_lloyd_taylor(T, E0, rb) - ER
    return residuals
def update_mergeseries(l6_info, output):
    """ Update the MergeSeries entry in l6_info."""
    source = l6_info["ERUsingLloydTaylor"]["outputs"][output]["source"]
    l6_info["Summary"]["EcosystemRespiration"].append(source)
    merge = l6_info["EcosystemRespiration"][source]["MergeSeries"]["source"].split(",")
    l6_info["MergeSeries"]["standard"][source] = {"output": source, "source": merge}
    return
def write_results(results, l6_info):
    """ Write the Lloyd-Taylor results to an Excel workbook."""
    file_path = l6_info["Files"]["file_path"]
    out_filename = l6_info["Files"]["out_filename"].replace(".nc", "_LloydTaylor.xlsx")
    xlsx_uri = os.path.join(file_path, out_filename)
    workbook = xlsxwriter.Workbook(xlsx_uri)
    write_results_settings(workbook, results, l6_info)
    write_results_summary(workbook, results)
    write_results_e0_windowed(workbook, results)
    write_results_rb_windowed(workbook, results)
    workbook.close()
    return
def write_results_settings(workbook, results, l6_info):
    """ Write the settings to the worksheet 'Settings'."""
    row = 0
    col = 0
    worksheet = workbook.add_worksheet("Settings")
    worksheet.write(row, col, "Parameters")
    outputs = list(l6_info["ERUsingLloydTaylor"]["outputs"].keys())
    items = ["target", "drivers", "minimum_temperature_spread", "minimum_window_points",
             "e0_prior", "rb_prior", "e0_step_size_days", "e0_window_size_days",
             "rb_step_size_days", "rb_window_size_days", "output_plots", "source"]
    for item in items:
        row = row + 1
        worksheet.write(row, col, item)
        for output in outputs:
            worksheet.write(row, col+1, output)
            if item in l6_info["ERUsingLloydTaylor"]["outputs"]["ER_LT_all"]:
                value = l6_info["ERUsingLloydTaylor"]["outputs"]["ER_LT_all"][item]
                if isinstance(value, list):
                    value = ",".join(value)
                worksheet.write(row, col+1, value)
    return
def write_results_summary(workbook, results):
    """ Write the annual E0 and rb values to the 'Summary' worksheet."""
    row = 0
    col = 0
    worksheet = workbook.add_worksheet("Summary")
    worksheet.write(row, col+1, "E0 year")
    worksheet.write(row, col+3, "rb year")
    worksheet.write(row, col+5, "E0 final")
    row = row + 1
    col = 0
    worksheet.write(row, col, "Year")
    for i in range(3):
        worksheet.write(row, i*2+1, "Value")
        worksheet.write(row, i*2+2, "Stderr")
    row = row + 1
    col = 0
    years = list(results["years"].keys())
    for year in years:
        res = results["years"][year]["all"]
        worksheet.write(row, col, year)
        worksheet.write(row, col+1, res["E0"]["value"])
        worksheet.write(row, col+2, res["E0"]["stderr"])
        worksheet.write(row, col+3, res["rb"]["value"])
        worksheet.write(row, col+4, res["rb"]["stderr"])
        rwqca = results["years"][year]["winqc_avg"]
        worksheet.write(row, col+5, rwqca["E0"]["value"])
        worksheet.write(row, col+6, rwqca["E0"]["stderr"])
        row = row + 1
    return
def write_results_e0_windowed(workbook, results):
    """ Write the E0 window results."""
    date_format = workbook.add_format({'num_format': 'yyyy-mm-dd HH:MM'})
    row = 0
    col = 0
    worksheet = workbook.add_worksheet("E0 windowed")
    worksheet.write(row, col+5, "E0")
    worksheet.write(row, col+7, "rb")
    row = row + 1
    worksheet.write(row, col, "Start")
    worksheet.write(row, col+1, "Mid")
    worksheet.write(row, col+2, "End")
    worksheet.write(row, col+3, "Window")
    worksheet.write(row, col+4, "Trange")
    worksheet.write(row, col+5, "Value")
    worksheet.write(row, col+6, "Stderr")
    worksheet.write(row, col+7, "Value")
    worksheet.write(row, col+8, "Stderr")
    row = row + 1
    years = list(results["years"].keys())
    for year in years:
        rw = results["years"][year]["win"]
        for k in range(len(rw["start"])):
            worksheet.write(row, col, rw["start"][k], date_format)
            worksheet.write(row, col+1, rw["mid"][k], date_format)
            worksheet.write(row, col+2, rw["end"][k], date_format)
            worksheet.write(row, col+3, rw["num"][k])
            worksheet.write(row, col+4, rw["Trange"][k])
            worksheet.write(row, col+5, rw["E0"]["value"][k])
            worksheet.write(row, col+6, rw["E0"]["stderr"][k])
            worksheet.write(row, col+7, rw["rb"]["value"][k])
            worksheet.write(row, col+8, rw["rb"]["stderr"][k])
            row = row + 1
    return
def write_results_rb_windowed(workbook, results):
    # write the rb window results
    date_format = workbook.add_format({'num_format': 'yyyy-mm-dd HH:MM'})
    row = 0
    col = 0
    worksheet = workbook.add_worksheet("Rb windowed")
    worksheet.write(row, col+5, "E0")
    worksheet.write(row, col+7, "rb")
    row = row + 1
    worksheet.write(row, col, "Start")
    worksheet.write(row, col+1, "Mid")
    worksheet.write(row, col+2, "End")
    worksheet.write(row, col+3, "Window")
    worksheet.write(row, col+4, "Trange")
    worksheet.write(row, col+5, "Value")
    worksheet.write(row, col+6, "Stderr")
    worksheet.write(row, col+7, "Value")
    worksheet.write(row, col+8, "Stderr")
    row = row + 1
    years = list(results["years"].keys())
    for year in years:
        rf = results["years"][year]["final"]
        for k in range(len(rf["start"])):
            worksheet.write(row, col, rf["start"][k], date_format)
            worksheet.write(row, col+1, rf["mid"][k], date_format)
            worksheet.write(row, col+2, rf["end"][k], date_format)
            worksheet.write(row, col+3, rf["num"][k])
            worksheet.write(row, col+4, rf["Trange"][k])
            worksheet.write(row, col+5, rf["E0"]["value"][k])
            worksheet.write(row, col+6, rf["E0"]["stderr"][k])
            worksheet.write(row, col+7, rf["rb"]["value"][k])
            worksheet.write(row, col+8, rf["rb"]["stderr"][k])
            row = row + 1
    return