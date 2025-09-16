# standard modules
import logging
import os
# 3rd party module
from configobj import ConfigObj
import numpy
# PFP modules
from scripts import constants as c
from scripts import pfp_io
from scripts import pfp_utils

pfp_log = os.environ["pfp_log"]
logger = logging.getLogger(pfp_log)

def generate_controlfile(main_ui):
    filename = pfp_io.get_filename_dialog(file_path=".", title='Choose a netCDF file')
    ds = pfp_io.NetCDFRead(filename)
    level = ds.root["Attributes"]["processing_level"]
    if level == "L1":
        generate_controlfile_l2(main_ui, ds)
    elif level == "L2":
        generate_controlfile_l3(main_ui, ds)
    else:
        msg = "Only L1 or L2 netCDF files can be used at present"
        logger.error(msg)
        return
    return
def generate_controlfile_l2(main_ui, ds):
    nc_url = ds.info["filepath"]
    l2_template_url = os.path.join("controlfiles", "templates", "L2", "L2_autogen_template.txt")
    cfg_template = ConfigObj(l2_template_url)
    cfg_template_labels = list(cfg_template["Variables"])
    av_labels = sorted(list(set([l.split("_")[0] for l in cfg_template_labels if "_Sd" not in l and "_Vr" not in l])))
    sd_labels = sorted(list(set([l.split("_")[0] for l in cfg_template_labels if "_Sd" in l])))
    vr_labels = sorted(list(set([l.split("_")[0] for l in cfg_template_labels if "_Vr" in l])))
    cv_labels = ["UxA", "UxC", "UxT", "UyA", "UyC", "UyT",
                 "UzA", "UzC", "UzT", "UxUy", "UxUz", "UyUz"]
    diag_labels = sorted(list(set([l for l in cfg_template_labels if "Diag" in l])))
    l2_template_labels = av_labels + sd_labels + vr_labels + cv_labels + diag_labels
    ds_labels = sorted(list(ds.root["Variables"]))
    if "DateTime" in ds_labels:
        ds_labels.remove("DateTime")
    cfg_site = ConfigObj(indent_type="    ", list_values=False)
    cfg_site["level"] = "L2"
    cfg_site["Files"] = {}
    cfg_site["Options"] = {}
    cfg_site["Variables"] = {}
    cfg_site["Plots"] = {}
    cfg_site["Files"]["file_path"] = os.path.dirname(nc_url)
    cfg_site["Files"]["in_filename"] = os.path.basename(nc_url)
    cfg_site["Files"]["out_filename"] = cfg_site["Files"]["in_filename"].replace("L1", "L2")
    cfg_site["Files"]["plot_path"] = os.path.join(cfg_site["Files"]["file_path"], "plots")
    cfg_site["Options"]["IRGA_Check"] = "Yes"
    cfg_site["Options"]["SONIC_Check"] = "Yes"
    # look for complete matches first
    for ds_label in list(ds_labels):
        if ds_label in cfg_template_labels:
            cfg_site["Variables"][ds_label] = cfg_template["Variables"][ds_label]
            ds_labels.remove(ds_label)
    for ds_label in list(ds_labels):
        ds_label_prefix = ds_label.split("_")[0]
        #if ds_label in av_labels or ds_label in sd_labels or ds_label in vr_labels or ds_label in diag_labels:
        if ds_label in l2_template_labels:
            if ds_label_prefix not in cfg_template_labels:
                continue
            cfg_site["Variables"][ds_label] = cfg_template["Variables"][ds_label_prefix]
            continue
        if ds_label_prefix in av_labels and ds_label_prefix in cfg_template_labels:
            cfg_site["Variables"][ds_label] = cfg_template["Variables"][ds_label_prefix]
    # adjust the range for presure based on the site altitude
    cfg_labels = list(cfg_site["Variables"].keys())
    if "altitude" in ds.root["Attributes"]:
        alt = str(ds.root["Attributes"]["altitude"])
        altitude = float(pfp_utils.strip_non_numeric(alt))
        mean_ps = 101.3 - 1.25*(altitude/100)
        upper_ps = mean_ps + 10
        lower_ps = mean_ps - 10
        ps_labels = [l for l in cfg_labels if l.split("_")[0] == "ps"]
        for ps_label in ps_labels:
            cfg_site["Variables"][ps_label]["RangeCheck"]["upper"] = str(upper_ps)
            cfg_site["Variables"][ps_label]["RangeCheck"]["lower"] = str(lower_ps)
    # set up the plot variables
    radiation_labels = ["Fsd", "Fsu", "Fld", "Flu", "Fn"]
    flux_labels = ["Fco2", "Fe", "Fh", "Fm", "ustar"]
    meteorology_labels = ["AH", "Ta", "RH", "Ws", "Wd", "ps", "Precip"]
    cov_sonic_labels = ["UxT", "UyT", "UzT", "UxUy", "UxUz", "UyUz"]
    cov_sonic_irga_labels = ["UxA", "UxC", "UyA", "UyC", "UzA", "UzC"]
    soil_temperature_labels = ["Ts"]
    soil_moisture_labels = ["Sws"]
    soil_heat_flux_labels = ["Fg"]
    # write the [Plots] section
    ds_labels = sorted(list(ds.root["Variables"].keys()))
    if "DateTime" in ds_labels:
        ds_labels.remove("DateTime")

    ds_labels_radiation = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in radiation_labels:
            ds_labels_radiation.append(ds_label)
    if len(ds_labels_radiation) > 0:
        cfg_site["Plots"]["Radiation"] = {"variables": ",".join(ds_labels_radiation)}

    ds_labels_cov_sonic = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in cov_sonic_labels:
            ds_labels_cov_sonic.append(ds_label)
    if len(ds_labels_cov_sonic) > 0:
        cfg_site["Plots"]["Covariances (SONIC)"] = {"variables": ",".join(ds_labels_cov_sonic)}

    ds_labels_cov_sonic_irga = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in cov_sonic_irga_labels:
            ds_labels_cov_sonic_irga.append(ds_label)
    if len(ds_labels_cov_sonic_irga) > 0:
        cfg_site["Plots"]["Covariances (SONIC & IRGA)"] = {"variables": ",".join(ds_labels_cov_sonic_irga)}

    ds_labels_flux = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in flux_labels:
            ds_labels_flux.append(ds_label)
    cfg_site["Plots"]["Turbulent fluxes"] = {"variables": ",".join(ds_labels_flux)}

    ds_labels_meteorology = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in meteorology_labels:
            ds_labels_meteorology.append(ds_label)
    if len(ds_labels_meteorology) > 0:
        cfg_site["Plots"]["Meteorology"] = {"variables": ",".join(ds_labels_meteorology)}

    ds_labels_soil_temperature = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_temperature_labels:
            ds_labels_soil_temperature.append(ds_label)
    st_labels = []
    for i in range(len(ds_labels_soil_temperature)):
        st = ds_labels_soil_temperature[i]
        if "cm" in st.split("_")[-1]:
            st_labels.append(st[0:st.index("cm")+2])
        elif "m" in st.split("_")[-1]:
            st_labels.append(st[0:st.index("m")+1])
        else:
            msg = "Unrecognised format for soil temperature variable name: " + st
            print(msg)
    st_labels = sorted(list(set(st_labels)))
    for n, st_label in enumerate(st_labels):
        plt_labels_soil_temperature = []
        for ds_label_soil_temperature in ds_labels_soil_temperature:
            if st_label in ds_label_soil_temperature:
                plt_labels_soil_temperature.append(ds_label_soil_temperature)
        if len(plt_labels_soil_temperature) > 0:
            cfg_site["Plots"]["Soil temperature: "+st_label] = {"variables": ",".join(plt_labels_soil_temperature)}

    ds_labels_soil_moisture = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_moisture_labels:
            ds_labels_soil_moisture.append(ds_label)
    sm_labels = []
    for i in range(len(ds_labels_soil_moisture)):
        sm = ds_labels_soil_moisture[i]
        if "cm" in sm.split("_")[-1]:
            sm_labels.append(sm[0:sm.index("cm")+2])
        elif "m" in sm.split("_")[-1]:
            sm_labels.append(sm[0:sm.index("m")+1])
        else:
            msg = "Unrecognised format for soil temperature variable name: " + sm
            print(msg)
    sm_labels = sorted(list(set(sm_labels)))
    for n, sm_label in enumerate(sm_labels):
        plt_labels_soil_moisture = []
        for ds_label_soil_moisture in ds_labels_soil_moisture:
            if sm_label in ds_label_soil_moisture:
                plt_labels_soil_moisture.append(ds_label_soil_moisture)
        #print(sm_label, len(plt_labels_soil_moisture), plt_labels_soil_moisture)
        if len(plt_labels_soil_moisture) > 0:
            cfg_site["Plots"]["Soil moisture: "+sm_label] = {"variables": ",".join(plt_labels_soil_moisture)}

    ds_labels_soil_heat_flux = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_heat_flux_labels:
            ds_labels_soil_heat_flux.append(ds_label)
    if len(ds_labels_soil_heat_flux) > 0:
        cfg_site["Plots"]["Soil heat flux"] = {"variables": ",".join(ds_labels_soil_heat_flux)}
    # save the automatically generated L2 control file
    cfg_filename = pfp_io.get_output_filename_dialog(file_path="L2_autogen.txt",
                                                     title="Choose an output file name ...")
    if len(str(cfg_filename)) == 0:
        return
    # write the control file
    cfg_site.filename = cfg_filename
    cfg_site.write()
    # open the control file in the PFP GUI
    main_ui.file_open(file_uri=cfg_filename)
    return
def generate_controlfile_l3(main_ui, ds):
    nc_url = ds.info["filepath"]
    ds_labels = list(ds.root["Variables"].keys())
    ds_gattr = ds.root["Attributes"].copy()
    if "canopy_height" in ds_gattr:
        canopy_height = pfp_utils.strip_non_numeric(ds_gattr["canopy_height"])
        canopy_height = round(float(canopy_height), 2)
    else:
        msg = " Global attribute 'canopy_height' missing, using default of 0.1 m"
        logger.warning(msg)
        canopy_height = 0.1
    got_height = False
    for label in ["CO2_IRGA_Av", "AH_IRGA_Av", "H2O_IRGA_Av"]:
        if label in ds_labels:
            var = pfp_utils.GetVariable(ds, label)
            if "height" in var["Attr"]:
                irga_height = pfp_utils.strip_non_numeric(var["Attr"]["height"])
                irga_height = round(float(irga_height), 2)
                got_height = True
                break
    if not got_height:
        msg = " No 'height' attribute found in CO2, AH or H2O variables, using default of 5 m"
        logger.warning(msg)
        irga_height = 5.0
    if "sonic_irga_separation_north" in ds_gattr:
        sonic_irga_separation_north = ds_gattr["sonic_irga_separation_north"]
        sonic_irga_separation_north = pfp_utils.strip_non_numeric(sonic_irga_separation_north)
        sonic_irga_separation_north = round(float(sonic_irga_separation_north), 3)
    else:
        msg = " Global attribute 'sonic_irga_separation_north' not found, using default of 0.1 m"
        logger.warning(msg)
        sonic_irga_separation_north = 0.1
    if "sonic_irga_separation_east" in ds_gattr:
        sonic_irga_separation_east = ds_gattr["sonic_irga_separation_east"]
        sonic_irga_separation_east = pfp_utils.strip_non_numeric(sonic_irga_separation_east)
        sonic_irga_separation_east = round(float(sonic_irga_separation_east), 3)
    else:
        msg = " Global attribute 'sonic_irga_separation_east' not found, using default of 0.1 m"
        logger.warning(msg)
        sonic_irga_separation_east = 0.1
    if "bulk_density" in ds_gattr:
        bulk_density = pfp_utils.strip_non_numeric(ds_gattr["bulk_density"])
        bulk_density = round(float(bulk_density), 0)
    else:
        msg = " Global attribute 'bulk_density' missing, using default of 1000 kg/m^3"
        logger.warning(msg)
        bulk_density = 1000
    if "organic_content" in ds_gattr:
        organic_content = pfp_utils.strip_non_numeric(ds_gattr["organic_content"])
        organic_content = round(float(organic_content), 2)
    else:
        msg = " Global attribute 'organic_content' missing, using default of 0.01"
        logger.warning(msg)
        organic_content = 0.01
    # create the control file
    cfg_site = ConfigObj(indent_type="    ", list_values=False)
    cfg_site["level"] = "L3"
    cfg_site["Files"] = {}
    cfg_site["Options"] = {}
    cfg_site["Soil"] = {}
    cfg_site["Massman"] = {}
    cfg_site["Variables"] = {}
    cfg_site["Plots"] = {}
    # populate the [Files] section
    cfg_site["Files"]["file_path"] = os.path.split(nc_url)[0]
    cfg_site["Files"]["in_filename"] = os.path.split(nc_url)[1]
    cfg_site["Files"]["out_filename"] = cfg_site["Files"]["in_filename"].replace("L2", "L3")
    cfg_site["Files"]["plot_path"] = os.path.join(cfg_site["Files"]["file_path"], "plots")
    # populate the [Options] section
    cfg_site["Options"]["ApplyFco2Storage"] = "No"
    cfg_site["Options"]["CorrectFgForStorage"] = "Yes"
    cfg_site["Options"]["CO2Units"] = "umol/mol"
    cfg_site["Options"]["Fco2Units"] = "umol/m^2/s"
    cfg_site["Options"]["KeepIntermediateSeries"] = "No"
    # disable calculation of fluxes from covariances unless all are present
    fluxes_need = ["Ux_SONIC_Av", "Uy_SONIC_Av", "Uz_SONIC_Av",
                   "UxUz", "UyUz", "UxUy",
                   "UzC", "UzA", "UzT",
                   "UxC", "UyC", "UxA",
                   "UyA", "UxT", "UyT"]
    got_needed = True
    for needed in fluxes_need:
        if needed not in ds_labels:
            got_needed = False
    fluxes_need_Sd = ["Ux_SONIC_Sd", "Uy_SONIC_Sd", "Uz_SONIC_Sd"]
    got_needed_Sd = True
    for needed in fluxes_need_Sd:
        if needed not in ds_labels:
            got_needed_Sd = False
    fluxes_need_Vr = ["Ux_SONIC_Vr", "Uy_SONIC_Vr", "Uz_SONIC_Vr"]
    got_needed_Vr = True
    for needed in fluxes_need_Vr:
        if needed not in ds_labels:
            got_needed_Vr = False
    if ((got_needed and got_needed_Sd) or (got_needed and got_needed_Vr)):
        cfg_site["Options"]["CalculateFluxes"] = "Yes"
        cfg_site["Options"]["MassmanCorrection"] = "Yes"
        cfg_site["Massman"]["north_separation"] = sonic_irga_separation_north
        cfg_site["Massman"]["east_separation"] = sonic_irga_separation_east
        cfg_site["Massman"]["zmd"] = round(irga_height - (0.66*canopy_height), 2)
    else:
        cfg_site["Options"]["CalculateFluxes"] = "No"
        del cfg_site["Massman"]
    # soil quantities will need to be in L2 global attributes
    cfg_site["Soil"]["BulkDensity"] = 1000
    cfg_site["Soil"]["OrganicContent"] = 0.01
    # meteorological data
    combine = "MergeSeries"
    tolerance = 0.1
    met_prefixes = ["AH", "CO2", "H2O", "RH", "Ta", "Wd", "Ws"]
    combine_close(cfg_site, ds, met_prefixes, combine=combine, tolerance=tolerance)
    cfgd_labels = cfg_site["Variables"].keys()
    if "Ta" in cfgd_labels:
        got_humidity = False
        for item in ["AH", "RH", "SH"]:
            if item in cfgd_labels:
                got_humidity = True
                break
        got_data = True
        for item in ["Tv_SONIC_Av", "ps"]:
            if item not in ds_labels:
                got_data = False
                break
        if got_humidity and got_data:
            src = cfg_site["Variables"]["Ta"][combine]["source"]
            cfg_site["Variables"]["Ta"][combine]["source"] = src + ",Ta_SONIC_Av"
    # radiation data
    combine = "MergeSeries"
    tolerance = 0.1
    rad_prefixes = ["Fn", "Fsd", "Fsu", "Fld", "Fld"]
    combine_close(cfg_site, ds, rad_prefixes, combine=combine, tolerance=tolerance)
    if "Fsd" in ds_labels and "Fsu" in ds_labels and "Fld" in ds_labels and "Flu" in ds_labels:
        if "Fn" in cfg_site["Variables"]:
            src = cfg_site["Variables"]["Fn"][combine]["source"]
            cfg_site["Variables"]["Fn"][combine]["source"] = "Fn_4cmpt," + src
        else:
            cfg_site["Variables"]["Fn"] = {combine: {"source": "Fn_4cmpt"}}
    # soil data
    combine = "AverageSeries"
    tolerance = 0.1
    soil_prefixes = ["Fg", "Sws", "Ts"]
    combine_close(cfg_site, ds, soil_prefixes, combine=combine, tolerance=tolerance)
    # flux data
    flux_prefixes = ["Fco2", "Fe", "Fh", "Fm", "ustar"]
    plt_flux_labels = []
    for prefix in flux_prefixes:
        if prefix not in cfg_site["Variables"]:
            cfg_site["Variables"][prefix] = {}
        prefix_labels = [l for l in ds_labels
                         if l.split("_")[0] == prefix
                         and l.split("_")[-1] in ["EP","EF"]]
        if len(prefix_labels) == 0:
            if cfg_site["Options"]["CalculateFluxes"] == "Yes":
                prefix_labels.append(prefix+"_PFP")
        else:
            cfg_site["Options"]["CalculateFluxes"] = "No"
            cfg_site["Options"]["MassmanCorrection"] = "No"
        cfg_site["Variables"][prefix]["MergeSeries"] = {"source": ",".join(prefix_labels)}
        plt_flux_labels.append(prefix)
    # remove variables that only have a MergeSeries key and the MergeSeries source is just the variable name
    cfg_site_labels = list(cfg_site["Variables"].keys())
    for cfg_site_label in cfg_site_labels:
        if "MergeSeries" in cfg_site["Variables"][cfg_site_label]:
            combine = "MergeSeries"
        elif "AverageSeries" in cfg_site["Variables"][cfg_site_label]:
            combine = "AverageSeries"
        else:
            continue
        src = cfg_site["Variables"][cfg_site_label][combine]["source"]
        if src == cfg_site_label:
            del cfg_site["Variables"][cfg_site_label][combine]
            if len(cfg_site["Variables"][cfg_site_label].keys()) == 0:
                del cfg_site["Variables"][cfg_site_label]
    # plots section
    radiation_labels = ["Fsd", "Fsu", "Fld", "Flu", "Fn"]
    flux_labels = ["Fco2", "Fe", "Fh", "Fm", "ustar"]
    meteorology_labels = ["AH", "Ta", "RH", "Ws", "Wd", "ps", "Precip"]
    soil_temperature_labels = ["Ts"]
    soil_moisture_labels = ["Sws"]
    soil_heat_flux_labels = ["Fg"]

    ds_labels = sorted(list(ds.root["Variables"].keys()))
    if "DateTime" in ds_labels:
        ds_labels.remove("DateTime")

    ds_labels_radiation = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in radiation_labels:
            ds_labels_radiation.append(ds_label)
    cfg_site["Plots"]["Radiation"] = {"variables": ",".join(ds_labels_radiation)}

    cfg_site["Plots"]["Turbulent fluxes"] = {"variables": ",".join(flux_labels)}

    cfg_site_labels = cfg_site["Variables"].keys()
    ds_labels_meteorology = []
    for label in meteorology_labels:
        if label in cfg_site_labels:
            ds_labels_meteorology.append(label)
            labels = [l for l in ds_labels if (l.split("_")[0] == label and
                                               l.split("_")[-1] not in ["Sd", "Vr"])]
            if len(labels) > 0:
                ds_labels_meteorology += labels
    cfg_site["Plots"]["Meteorology"] = {"variables": ",".join(ds_labels_meteorology)}

    ds_labels_soil_temperature = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_temperature_labels:
            ds_labels_soil_temperature.append(ds_label)
    st_labels = []
    for i in range(len(ds_labels_soil_temperature)):
        st = ds_labels_soil_temperature[i]
        if "cm" in st.split("_")[-1]:
            st_labels.append(st[0:st.index("cm")+2])
        elif "m" in st.split("_")[-1]:
            st_labels.append(st[0:st.index("m")+1])
        else:
            msg = "Unrecognised format for soil temperature variable name: " + st
            print(msg)
    st_labels = sorted(list(set(st_labels)))
    for n, st_label in enumerate(st_labels):
        plt_labels_soil_temperature = [st_label]
        for ds_label_soil_temperature in ds_labels_soil_temperature:
            if st_label in ds_label_soil_temperature:
                plt_labels_soil_temperature.append(ds_label_soil_temperature)
        cfg_site["Plots"]["Soil temperature: "+st_label] = {"variables": ",".join(plt_labels_soil_temperature)}

    ds_labels_soil_moisture = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_moisture_labels:
            ds_labels_soil_moisture.append(ds_label)
    sm_labels = []
    for i in range(len(ds_labels_soil_moisture)):
        sm = ds_labels_soil_moisture[i]
        if "cm" in sm.split("_")[-1]:
            sm_labels.append(sm[0:sm.index("cm")+2])
        elif "m" in sm.split("_")[-1]:
            sm_labels.append(sm[0:sm.index("m")+1])
        else:
            msg = "Unrecognised format for soil temperature variable name: " + sm
            print(msg)
    sm_labels = sorted(list(set(sm_labels)))
    for n, sm_label in enumerate(sm_labels):
        plt_labels_soil_moisture = [sm_label]
        for ds_label_soil_moisture in ds_labels_soil_moisture:
            if sm_label in ds_label_soil_moisture:
                plt_labels_soil_moisture.append(ds_label_soil_moisture)
        cfg_site["Plots"]["Soil moisture: "+sm_label] = {"variables": ",".join(plt_labels_soil_moisture)}

    ds_labels_soil_heat_flux = []
    for ds_label in ds_labels:
        ds_label_prefix = ds_label.split("_")[0]
        if ds_label_prefix in soil_heat_flux_labels:
            if ds_label_prefix not in ds_labels_soil_heat_flux:
                ds_labels_soil_heat_flux.append(ds_label_prefix)
            ds_labels_soil_heat_flux.append(ds_label)
    cfg_site["Plots"]["Soil heat flux"] = {"variables": ",".join(ds_labels_soil_heat_flux)}
    # save the automatically generated L2 control file
    cfg_filename = pfp_io.get_output_filename_dialog(file_path="L3_autogen.txt",
                                                     title="Choose an output file name ...")
    if len(str(cfg_filename)) == 0:
        return
    # write the control file
    cfg_site.filename = cfg_filename
    cfg_site.write()
    # open the control file in the PFP GUI
    main_ui.file_open(file_uri=cfg_filename)
    return
def substrings(s):
    return {s[i:j]
            for j in range(len(s)+1)
            for i in range(j+1)}
def combine_close(cfgd, ds, prefixes, combine="MergeSeries", tolerance=0.1):
    ds_labels = list(ds.root["Variables"].keys())
    for prefix in prefixes:
        prefix_labels = [l for l in ds_labels if l.split("_")[0] == prefix]
        if len(prefix_labels) == 0:
            continue
        # get dictionaries of variables and variable information
        td, variables = get_dicts(ds, prefix_labels)
        # get groups of variables that are close in height/depth
        get_close_variables(td, tolerance)
        # find highest measurements, create an entry for "prefix" and get the combine order
        # for 'fast' instruments
        get_combine_order(cfgd, variables, td, prefix, combine)
        # make the variable entries in the control file
        get_variable_entries(cfgd, td, combine)
    return
def get_dicts(ds, prefix_labels):
    td = {}
    variables = {}
    for prefix_label in prefix_labels:
        try:
            variables[prefix_label] = pfp_utils.GetVariable(ds, prefix_label)
            if variables[prefix_label]["Attr"]["statistic_type"] != "average":
                del variables[prefix_label]
                continue
            td[prefix_label] = {"height": 0, "instrument": "", "coverage": 0, "close": []}
            ha = variables[prefix_label]["Attr"]["height"]
            if "to" in ha:
                n = len(ha.split("to"))
                if n == 2:
                    hl = [float(pfp_utils.strip_non_numeric(s.strip())) for s in ha.split("to")]
                    h = round(sum(hl)/len(hl), 2)
                else:
                    msg = "Unrecognised 'height' attribute format (" + prefix_label + ")"
                    logger.warning(msg)
                    del td[prefix_label]
                    continue
            else:
                h = float(pfp_utils.strip_non_numeric(ha))
            td[prefix_label]["height"] = h
            td[prefix_label]["instrument"] = variables[prefix_label]["Attr"]["instrument"]
            coverage = round(100*numpy.ma.count(variables[prefix_label]["Data"])/len(variables[prefix_label]["Data"]))
            td[prefix_label]["coverage"] = coverage
        except Exception as e:
            print(prefix_label + " gave the following exception")
            print(type(e), e)
            del td[prefix_label]
            continue
    return td, variables
def get_close_variables(td, tolerance):
    used = []
    # find all variables that are within 'tolerance' of the same height
    all_keys = list(td.keys())
    for key1 in all_keys:
        if key1 in used:
            continue
        height1 = td[key1]["height"]
        for key2 in all_keys:
            height2 = td[key2]["height"]
            if abs(height2 - height1) < tolerance*abs(height1):
                td[key1]["close"].append(key2)
                used.append(key2)
    # delete empty entries
    all_keys = list(td.keys())
    for key in all_keys:
        if key in td:
            if len(td[key]["close"]) == 0:
                del td[key]
    return
def get_combine_order(cfgd, variables, td, prefix, combine):
    open_path_irgas = list(c.instruments["irgas"]["open_path"].keys())
    closed_path_irgas = list(c.instruments["irgas"]["closed_path"].keys())
    irgas = open_path_irgas + closed_path_irgas
    sonics = list(c.instruments["sonics"].keys())
    keys = list(td.keys())
    heights = []
    for key in keys:
        heights.append(td[key]["height"])
    max_height_index = heights.index(max(heights))
    label_highest = keys[max_height_index]
    td[label_highest][combine] = td[label_highest]["close"].copy()
    instruments = [variables[l]["Attr"]["instrument"] for l in td[label_highest]["close"]]
    last_labels = []
    for label, instrument in zip(td[label_highest][combine], instruments):
        if instrument in irgas+sonics:
            last_labels.append(label)
    for last_label in last_labels:
        td[label_highest][combine].remove(last_label)
        td[label_highest][combine].append(last_label)
    td[prefix] = td[label_highest].copy()
    if prefix not in cfgd["Variables"]:
        cfgd["Variables"][prefix] = {}
    cfgd["Variables"][prefix][combine] = {"source": ",".join(td[prefix][combine])}
    del td[prefix]
    return
def get_variable_entries(cfgd, td, combine):
    for item in list(td.keys()):
        if len(td[item]["close"]) > 1:
            common = set.intersection(*map(substrings, td[item]["close"]))
            label = max(common, key=len)
            if label[-1] != "m":
                continue
        elif len(td[item]["close"]) == 1:
            l = td[item]["close"][0]
            p = l.split("_")
            if "cm" in p[-1]:
                label = l[:l.index("cm")+2]
            elif "m" in p[-1]:
                label = l[:l.index("m")+1]
            else:
                label = l
        if label not in cfgd["Variables"]:
            cfgd["Variables"][label] = {}
        cfgd["Variables"][label][combine] = {"source": ",".join(td[item]["close"])}
    return