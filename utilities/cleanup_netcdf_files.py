# standard modules
import datetime
import copy
import logging
import os
import subprocess
import sys
# 3rd party modules
from configobj import ConfigObj
import numpy
from PyQt5 import QtWidgets
import timezonefinder
import xlrd
# PFP modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("utilities")-1])
import scripts.constants as c
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_ts as pfp_ts
import scripts.pfp_utils as pfp_utils

now = datetime.datetime.now()
log_file_name = "cleanup_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_name = os.path.join("logfiles", log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_name, to_screen=True)
cfchecker_file_name = log_file_name.replace(".log", "_cfchecker.all")
cfchecker_file = open(cfchecker_file_name, "w")
error_file_name = cfchecker_file_name.replace(".all", ".errors")
warning_file_name = cfchecker_file_name.replace(".all", ".warnings")

def change_global_attributes(std, ds):
    """
    Purpose:
     Clean up the global attributes.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # check site_name is in ds.root["Attributes"]
    gattr_list = list(ds.root["Attributes"].keys())
    if "site_name" not in gattr_list:
        msg = "Global attributes: site_name not found"
        logger.warning(msg)
    # check latitude and longitude are in ds.root["Attributes"]
    if "latitude" not in gattr_list:
        msg = "Global attributes: latitude not found"
        logger.warning(msg)
    else:
        lat_string = str(ds.root["Attributes"]["latitude"])
        if len(lat_string) == 0:
            msg = "Global attributes: latitude empty"
            logger.warning(msg)
        else:
            lat = pfp_utils.convert_anglestring(lat_string)
        ds.root["Attributes"]["latitude"] = str(lat)
    if "longitude" not in gattr_list:
        msg = "Global attributes: longitude not found"
        logger.warning(msg)
    else:
        lon_string = str(ds.root["Attributes"]["longitude"])
        if len(lon_string) == 0:
            msg = "Global attributes: longitude empty"
            logger.warning(msg)
        else:
            lon = pfp_utils.convert_anglestring(lon_string)
        ds.root["Attributes"]["longitude"] = str(lon)
    # check to see if there there is a time_zone global attribute
    gattr_list = list(ds.root["Attributes"].keys())
    if not "time_zone" in gattr_list:
        # get the site name
        site_name = ds.root["Attributes"]["site_name"]
        sn = site_name.replace(" ","").replace(",","").lower()
        # first, see if the site is in constants.tz_dict
        if sn in list(c.tz_dict.keys()):
            ds.root["Attributes"]["time_zone"] = c.tz_dict[sn]
        else:
            if "latitude" in gattr_list and "longitude" in gattr_list:
                lat = float(ds.root["Attributes"]["latitude"])
                lon = float(ds.root["Attributes"]["longitude"])
                if lat != -9999 and lon != -9999:
                    tf = timezonefinder.TimezoneFinder()
                    tz = tf.timezone_at(lng=lon, lat=lat)
                    ds.root["Attributes"]["time_zone"] = tz
                else:
                    msg = "Global attributes: unable to define time zone"
                    logger.warning(msg)
                    ds.root["Attributes"]["time_zone"] = ""
    # remove deprecated global attributes
    flag_list = [g for g in list(ds.root["Attributes"].keys()) if "Flag" in g]
    others_list = pfp_utils.string_to_list(std["Global_attributes"]["deprecated"]["global"])
    remove_list = others_list + flag_list
    for gattr in list(ds.root["Attributes"].keys()):
        if gattr in remove_list:
            ds.root["Attributes"].pop(gattr)
    # rename global attributes
    for item in std["Global_attributes"]["rename_exact"]:
        if item in list(ds.root["Attributes"].keys()):
            new_key = std["Global_attributes"]["rename_exact"][item]
            ds.root["Attributes"][new_key] = ds.root["Attributes"].pop(item)
    # replace space characters with underscore
    gattrs = sorted(list(ds.root["Attributes"].keys()))
    for gattr in gattrs:
        if " " in gattr:
            new_gattr = gattr.replace(" ", "_")
            ds.root["Attributes"][new_gattr] = ds.root["Attributes"].pop(gattr)
    # add or change global attributes as required
    gattr_list = sorted(list(std["Global_attributes"]["force"].keys()))
    for gattr in gattr_list:
        ds.root["Attributes"][gattr] = std["Global_attributes"]["force"][gattr]
    # add acknowledgement if not present
    gattrs = sorted(list(ds.root["Attributes"].keys()))
    if "acknowledgement" not in gattrs:
        if "acknowledgement" in std["Global_attributes"]:
            msg = std["Global_attributes"]["acknowledgement"]
        else:
            msg = "This work used eddy covariance data collected by the TERN Ecosystem "
            msg += "Processes facility. Ecosystem Processes would like to acknowledge the financial support of the "
            msg += "Australian Federal Government via the National Collaborative Research Infrastructure Scheme "
            msg += "and the Education Investment Fund."
        ds.root["Attributes"]["acknowledgement"] = msg
    # force the fluxnet_id
    site_name = ds.root["Attributes"]["site_name"].replace(" ","")
    if site_name in std["Global_attributes"]["fluxnet_id"]:
        ds.root["Attributes"]["fluxnet_id"] = std["Global_attributes"]["fluxnet_id"][site_name]
    return

def change_variable_attributes(std, ds):
    """
    Purpose:
     Clean up the variable attributes.
    Usage:
    Author: PRI
    Date: November 2018
    """
    # coerce units into a standard form
    msg = " Parse variable attributes"
    logger.info(msg)
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        # parse variable attributes to new format
        variable["Attr"] = parse_variable_attributes(variable["Attr"])
        pfp_utils.CreateVariable(ds, variable)

    # coerce units into a standard form
    msg = " Changing variable units"
    logger.info(msg)
    old_units = list(std["Variables"]["units_map"].keys())
    deprecated_units = pfp_utils.string_to_list(std["Variables"]["deprecated"]["units"])
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        old_unit = variable["Attr"]["units"]
        if old_unit in old_units:
            variable["Attr"]["units"] = std["Variables"]["units_map"][old_unit]
            pfp_utils.CreateVariable(ds, variable)
        elif old_unit in deprecated_units or len(old_unit) == 0:
            variable["Attr"]["units"] = "1"
            pfp_utils.CreateVariable(ds, variable)

    # do any units conversion required
    msg = " Converting units"
    logger.info(msg)
    stdv = std["Variables"]
    stdva = std["Variables"]["attributes"]
    labels = list(ds.root["Variables"].keys())
    vattrs = pfp_utils.string_to_list(stdv["units_convert"]["labels"])
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        for vattr in vattrs:
            #if label[:len(vattr)] == vattr:
            if label.split("_")[0] == vattr:
                old_units = variable["Attr"]["units"]
                std_units = stdva[vattr]["units"]
                if old_units != std_units:
                    msg = " Converting " + label + " to units of " + std_units
                    logger.info(msg)
                    variable = pfp_utils.convert_units_func(ds, variable, std_units)
                    pfp_utils.CreateVariable(ds, variable)

    # rename existing long_name to description, introduce a
    # consistent long_name attribute
    msg = " Long name, description and units"
    logger.info(msg)
    stdva = std["Variables"]["attributes"]
    vattr_list = list(stdva.keys())
    labels = list(ds.root["Variables"].keys())
    descr = "description_" + ds.root["Attributes"]["processing_level"]
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        # existing long_name to description
        variable["Attr"][descr] = variable["Attr"]["long_name"]
        for item in vattr_list:
            if label == item:
                for key in list(stdva[item].keys()):
                    if key == "units":
                        if "_Sd" not in label and "_Vr" not in label:
                            old_units = variable["Attr"]["units"]
                            new_units = pfp_utils.string_to_list(stdva[item]["units"])
                            if old_units not in new_units and old_units != "1":
                                msg = "Units '" + old_units + "' not found for variable " + label
                                logger.warning(msg)
                                continue
                            if old_units in new_units:
                                if "standard_name" in stdva[item]:
                                    idx = new_units.index(old_units)
                                    standard_names = pfp_utils.string_to_list(stdva[item]["standard_name"])
                                    variable["Attr"]["standard_name"] = standard_names[idx]
                    elif key == "standard_name":
                        if "_Sd" in label or "_Vr" in label:
                            # remove standard_name for standard deviations and variances
                            if "standard_name" in variable["Attr"]:
                                variable["Attr"].pop("standard_name")
                    else:
                        variable["Attr"][key] = stdva[item][key]
                # remove standard_name attribute if not defined for this variable
                if (("standard_name" in variable["Attr"]) and ("standard_name" not in stdva[item])):
                    variable["Attr"].pop("standard_name")
            elif label.split("_")[0] == item:
                for key in list(stdva[item].keys()):
                    if key == "units":
                        if "_Sd" not in label and "_Vr" not in label:
                            old_units = variable["Attr"]["units"]
                            new_units = pfp_utils.string_to_list(stdva[item]["units"])
                            if old_units not in new_units and old_units != "1":
                                msg = "Units '" + old_units + "' not found for variable " + label
                                logger.warning(msg)
                                continue
                            if old_units in new_units:
                                if "standard_name" in stdva[item]:
                                    idx = new_units.index(old_units)
                                    standard_names = pfp_utils.string_to_list(stdva[item]["standard_name"])
                                    variable["Attr"]["standard_name"] = standard_names[idx]
                    elif key == "standard_name":
                        if "_Sd" in label or "_Vr" in label:
                            # remove standard_name for standard deviations and variances
                            if "standard_name" in variable["Attr"]:
                                variable["Attr"].pop("standard_name")
                    else:
                        variable["Attr"][key] = stdva[item][key]
                # remove standard_name attribute if not defined for this variable
                if (("standard_name" in variable["Attr"]) and ("standard_name" not in stdva[item])):
                    variable["Attr"].pop("standard_name")
        pfp_utils.CreateVariable(ds, variable)

    # parse variable attributes to new format, remove deprecated variable attributes
    # and fix valid_range == "-1e+35,1e+35"
    msg = " Remove deprecated, check valid_range"
    logger.info(msg)
    tmp = std["Variables"]["deprecated"]["attributes"]
    deprecated = pfp_utils.string_to_list(tmp)
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        ## parse variable attributes to new format
        #variable["Attr"] = parse_variable_attributes(variable["Attr"])
        # remove deprecated variable attributes
        for vattr in deprecated:
            if vattr in list(variable["Attr"].keys()):
                del variable["Attr"][vattr]
        # fix valid_range == "-1e+35,1e+35"
        if "valid_range" in variable["Attr"]:
            valid_range = variable["Attr"]["valid_range"]
            if valid_range == "-1e+35,1e+35":
                d = numpy.ma.min(variable["Data"])
                mn = pfp_utils.round2significant(d, 4, direction='down')
                d = numpy.ma.max(variable["Data"])
                mx = pfp_utils.round2significant(d, 4, direction='up')
                variable["Attr"]["valid_range"] = repr(mn) + "," + repr(mx)
        pfp_utils.CreateVariable(ds, variable)

    # ugly hack to deal with CO2_IRGA_Av, CO2_IRGA_Sd and CO2_IRGA_Vr units
    msg = " Converting CO2 units"
    logger.info(msg)
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        if (label.split("_")[0] == "CO2"):
            co2 = pfp_utils.GetVariable(ds, label)
            if ((label == "CO2") and (co2["Attr"]["units"] == "mg/m^3")):
                co2 = pfp_utils.convert_units_func(ds, co2, "umol/mol")
                pfp_utils.CreateVariable(ds, co2)
            elif ((label[-3:] == "_Av") and (co2["Attr"]["units"] == "mg/m^3")):
                co2 = pfp_utils.convert_units_func(ds, co2, "umol/mol")
                pfp_utils.CreateVariable(ds, co2)
            elif ((label[-3:] == "_Sd") and (co2["Attr"]["units"] == "mg/m^3")):
                co2 = pfp_utils.convert_units_func(ds, co2, "mmol/m^3")
                # convert_units_func() will add a standard name when units=mmol/m^3
                # here we dump standard_name for a standard deviation
                if "standard_name" in co2["Attr"]:
                    co2["Attr"].pop("standard_name")
                pfp_utils.CreateVariable(ds, co2)
            elif ((label[-3:] == "_Vr") and (co2["Attr"]["units"] == "mg^2/m^6")):
                msg = " Converting " + label + " from " + co2["Attr"]["units"]
                msg += " to mmol^2/m^6"
                logger.info(msg)
                co2["Data"] = numpy.ma.sqrt(co2["Data"])
                co2["Attr"]["units"] = "mg/m^3"
                co2 = pfp_utils.convert_units_func(ds, co2, "mmol/m^3")
                co2["Data"] = co2["Data"]*co2["Data"]
                co2["Attr"]["units"] = "mmol^2/m^6"
                # convert_units_func() will add a standard name when units=mmol/m^3
                # here we dump standard_name for a variance
                if "standard_name" in co2["Attr"]:
                    co2["Attr"].pop("standard_name")
                pfp_utils.CreateVariable(ds, co2)
            else:
                if co2["Attr"]["units"] not in ["umol/mol", "mmol/m^3", "mmol^2/m^6", "1"]:
                    msg = " Unrecognised units for " + label
                    msg += " (" + co2["Attr"]["units"] + ")"
                    logger.warning(msg)

    # set the statistic_type variable attribute for standard deviations and
    # variances and remove standard_name if present
    msg = " Setting statistic_type"
    logger.info(msg)
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        if label[-3:] == "_Sd":
            variable["Attr"]["statistic_type"] = "standard_deviation"
            if "standard_name" in variable["Attr"]:
                variable["Attr"].pop("standard_name")
            pfp_utils.CreateVariable(ds, variable)
        elif label[-3:] == "_Vr":
            variable["Attr"]["statistic_type"] = "variance"
            if "standard_name" in variable["Attr"]:
                variable["Attr"].pop("standard_name")
            pfp_utils.CreateVariable(ds, variable)
        else:
            pass

    # remove undefined attributes
    msg = " Remove undefined attributes"
    logger.info(msg)
    stdvd = std["Variables"]["deprecated"]
    units_deprecated = pfp_utils.string_to_list(stdvd["units"])
    height_deprecated = pfp_utils.string_to_list(stdvd["height"])
    miscellaneous_deprecated = pfp_utils.string_to_list(stdvd["miscellaneous"])
    standard_name_deprecated = pfp_utils.string_to_list(stdvd["standard_name"])
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        attrs = list(variable["Attr"].keys())
        for attr in attrs:
            if attr == "units":
                if ((variable["Attr"][attr] in units_deprecated) or
                    (len(str(variable["Attr"][attr])) == 0)):
                    variable["Attr"].pop(attr)
            elif attr == "height":
                if ((variable["Attr"][attr] in height_deprecated) or
                    (len(str(variable["Attr"][attr])) == 0)):
                    variable["Attr"].pop(attr)
            elif attr == "standard_name":
                if ((variable["Attr"][attr] in standard_name_deprecated) or
                    (len(str(variable["Attr"][attr])) == 0)):
                    variable["Attr"].pop(attr)
            elif len(str(variable["Attr"][attr])) == 0:
                variable["Attr"].pop(attr)
            else:
                if ((variable["Attr"][attr] in miscellaneous_deprecated) or
                    (len(str(variable["Attr"][attr])) == 0)):
                    variable["Attr"].pop(attr)
        pfp_utils.CreateVariable(ds, variable)

    # append '%' symbol to coverage_L2 and coverage_L3 variable sttributes
    msg = " Appending % to coverage statistics"
    logger.info(msg)
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        variable = pfp_utils.GetVariable(ds, label)
        attrs = list(variable["Attr"].keys())
        for attr in attrs:
            if "coverage_" in attr:
                if "%" not in variable["Attr"][attr]:
                    variable["Attr"][attr] += "%"
        pfp_utils.CreateVariable(ds, variable)
    return

def change_variable_names(std, ds):
    """
    Purpose:
     Change variable names to the new (October 2018) scheme.
    Usage:
    Author: PRI
    Date: October 2018
    """
    msg = " Changing variable names"
    logger.info(msg)
    # get a list of exact mappings
    renames_exact = list(std["Variables"]["rename_exact"].keys())
    # loop over the variables in the data structure
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        if label in renames_exact:
            new_name = std["Variables"]["rename_exact"][label]
            ds.root["Variables"][new_name] = ds.root["Variables"].pop(label)
    # get a list of pattern mappings
    renames_pattern = list(std["Variables"]["rename_pattern"].keys())
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        for rp in renames_pattern:
            #if label[:len(rp)] == rp:
                #np = std["Variables"]["rename_pattern"][rp]
                #new_label = label.replace(rp, np)
                #ds.root["Variables"][new_label] = ds.root["Variables"].pop(label)
            if label.split("_")[0] == rp:
                np = std["Variables"]["rename_pattern"][rp]
                new_label = label.replace(rp, np)
                ds.root["Variables"][new_label] = ds.root["Variables"].pop(label)
    # replace "." in variable names with "p"
    labels = list(ds.root["Variables"].keys())
    for label in labels:
        if "." in label:
            new_label = label.replace(".", "p")
            ds.root["Variables"][new_label] = ds.root["Variables"].pop(label)
    return

def consistent_Fco2_storage(std, ds, site):
    """
    Purpose:
     Make the various incarnations of single point Fco2 storage consistent.
    Author: PRI
    Date: November 2019
    """
    msg = " Sorting Fco2_storage ..."
    logger.info(msg)
    ## save Fc_single if it exists - debug only
    #labels = ds.root["Variables"].keys()
    #if "Fc_single" in labels:
        #variable = pfp_utils.GetVariable(ds, "Fc_single")
        #variable["Label"] = "Fc_sinorg"
        #pfp_utils.CreateVariable(ds, variable)
        #pfp_utils.DeleteVariable(ds, "Fc_single")
    # do nothing if Fco2_single exists
    labels = list(ds.root["Variables"].keys())
    if "Fco2_single" in labels:
        pass
    # Fco2_single may be called Fco2_storage
    elif "Fco2_storage" in labels:
        level = ds.root["Attributes"]["processing_level"]
        descr = "description_" + level
        variable = pfp_utils.GetVariable(ds, "Fco2_storage")
        if (("single" in variable["Attr"][descr]) or
            (site in ["AdelaideRiver"])):
            variable["Label"] = "Fco2_single"
            pfp_utils.CreateVariable(ds, variable)
            pfp_utils.DeleteVariable(ds, "Fco2_storage")
    else:
        # neither Fco2_single nor Fco2_storage exist, try to calculate
        # check to see if the measurement height is defined
        zms = None
        CO2 = pfp_utils.GetVariable(ds, "CO2")
        if "height" in CO2["Attr"]:
            zms = pfp_utils.get_number_from_heightstring(CO2["Attr"]["height"])
        if zms is None:
            xls_name = std["Files"]["site_information"]
            site_information = xl_read_site_information(xls_name, site)
            if len(site_information) != 0:
                s = site_information["IRGA"]["Height"]
                zms = pfp_utils.get_number_from_heightstring(s)
            else:
                while zms is None:
                    file_name = std["Files"]["in_filename"]
                    prompt = "Enter CO2 measuement height in metres"
                    text, ok = QtWidgets.QInputDialog.getText(None, file_name,
                                                              prompt,
                                                              QtWidgets.QLineEdit.Normal,
                                                              "")
                    zms = pfp_utils.get_number_from_heightstring(text)
        # update the CO2 variable attribute
        CO2["Attr"]["height"] = zms
        pfp_utils.CreateVariable(ds, CO2)
        # calculate single point Fc storage term
        cf = {"Options": {"zms": zms}}
        info = {"CO2": {"label": "CO2", "height": zms}}
        pfp_ts.CalculateFco2StorageSinglePoint(cf, ds, info)
        # convert Fco2_single from mg/m2/s to umol/m2/s
        pfp_utils.CheckUnits(ds, "Fco2_single", "umol/m^2/s", convert_units=True)
    return

def copy_ws_wd(ds):
    """
    Purpose:
     Make sure the Ws and Wd variables are in the L3 netCDF files.
    Usage:
    Author: PRI
    Date: October 2018
    """
    msg = " Sorting Ws and Wd ..."
    logger.info(msg)
    # get a list of the series
    series_list = sorted(list(ds.root["Variables"].keys()))
    if "Wd" not in series_list:
        if "Wd_SONIC_Av" in series_list:
            ds.root["Variables"]["Wd"] = copy.deepcopy(ds.root["Variables"]["Wd_SONIC_Av"])
            ds.root["Variables"]["Wd"]["Attr"]["long_name"] = "Wind direction (copied from Wd_SONIC_Av)"
    if "Ws" not in series_list:
        if "Ws_SONIC_Av" in series_list:
            ds.root["Variables"]["Ws"] = copy.deepcopy(ds.root["Variables"]["Ws_SONIC_Av"])
            ds.root["Variables"]["Ws"]["Attr"]["long_name"] = "Wind speed (copied from Ws_SONIC_Av)"
    return

def exclude_variables(std, ds):
    """
    Purpose:
     Remove deprecated variables from a netCDF file.
    Usage:
    Author: PRI
    Date: October 2018
    """
    msg = " Excluding variables ..."
    logger.info(msg)
    series_list = sorted(list(ds.root["Variables"].keys()))
    exclude_list = pfp_utils.string_to_list(std["Variables"]["exclude"]["exclude"])
    flag_list = [v+"_QCFlag" for v in exclude_list if v+"_QCFlag" in series_list]
    remove_list = exclude_list + flag_list
    for label in series_list:
        if label in remove_list:
            ds.root["Variables"].pop(label)
    return

def include_variables(std, ds_in):
    """
    Purpose:
     Only pick variables that match the specified string for the length
     of the specified string.
    Usage:
    Author: PRI
    Date: November 2018
    """
    msg = " Including variables ..."
    logger.info(msg)
    # get a new data structure
    ds_out = pfp_io.DataStructure()
    # copy the global attributes
    for gattr in ds_in.root["Attributes"]:
        ds_out.root["Attributes"][gattr] = ds_in.root["Attributes"][gattr]
    # loop over variables to be included
    include_list = pfp_utils.string_to_list(std["Variables"]["include"]["include"])
    series_list = list(ds_in.root["Variables"].keys())
    for item in include_list:
        for label in series_list:
            if label[0:len(item)] == item:
                ds_out.root["Variables"][label] = ds_in.root["Variables"][label]
    return ds_out

def parse_variable_attributes(attributes):
    """
    Purpose:
     Clean up the variable attributes.
    Usage:
    Author: PRI
    Date: September 2019
    """
    for attr in attributes:
        value = attributes[attr]
        if not isinstance(value, str):
            continue
        if attr in ["rangecheck_lower", "rangecheck_upper", "diurnalcheck_numsd"]:
            if ("[" in value) and ("]" in value) and ("*" in value):
                # old style of [value]*12
                value = value[value.index("[")+1:value.index("]")]
            elif ("[" in value) and ("]" in value) and ("*" not in value):
                # old style of [1,2,3,4,5,6,7,8,9,10,11,12]
                value = value.replace("[", "").replace("]", "")
            strip_list = [" ", '"', "'"]
        elif ("ExcludeDates" in attr or
              "ExcludeHours" in attr or
              "LowerCheck" in attr or
              "UpperCheck" in attr):
            strip_list = ["[", "]", '"', "'"]
        else:
            strip_list = ['"', "'"]
        for c in strip_list:
            if c in value:
                value = value.replace(c, "")
        attributes[attr] = value
    return attributes

def xl_read_site_information(xls_name, site_name):
    """
    Purpose:
     Read the site information workbook.
    Usage:
     site_information = pfp_io.xl_read_site_information(xls_name)
    Author: PRI
    Date: December 2019
    """
    xl_book = xlrd.open_workbook(xls_name)
    sheet_names = [sn.replace(" ", "") for sn in xl_book.sheet_names()]
    if site_name.replace(" ", "") not in sheet_names:
        msg = " Site " + str(site_name.replace(" ", ""))
        msg += " not found in site information workbook"
        logger.warning(msg)
        return {}
    xl_sheet = xl_book.sheet_by_index(sheet_names.index(site_name.replace(" ", "")))
    nrows = xl_sheet.nrows
    ncols = xl_sheet.ncols
    info = {"site_name": site_name}
    for row in range(1, nrows):
        measurement = xl_sheet.cell_value(row, 0)
        info[measurement] = {}
        for col in range(1, ncols):
            field = xl_sheet.cell_value(0, col)
            info[measurement][field] = xl_sheet.cell_value(row, col)
    return info

logger = logging.getLogger("pfp_log")

std_name = os.path.join("..", "controlfiles", "standard", "update_control_files.txt")
if os.path.exists(std_name):
    std = ConfigObj(std_name, indent_type="    ", list_values=False, write_empty_values=True)
else:
    msg = " 'update_control_files.txt' control file not found"
    logger.error(msg)

rp = os.path.join(os.sep, "mnt", "OzFlux", "Sites")
#rp = os.path.join(os.sep, "home", "peter", "WD2TB", "OzFlux", "Sites")
#rp = os.path.join(os.sep, "home", "peter", "OzFlux", "Sites")
#sites = ["DalyUncleared"]
sites = ["AdelaideRiver", "AliceSpringsMulga", "Boyagin", "Calperum", "CapeTribulation", "Collie",
         "CowBay", "CumberlandPlain", "DalyPasture", "DalyUncleared", "DryRiver", "Emerald",
         "Fletcherview", "FoggDam", "Gingin", "GreatWesternWoodlands", "HowardSprings", "Litchfield",
         "Longreach", "Loxton", "Otway", "RedDirtMelonFarm", "Ridgefield", "RiggsCreek", "RobsonCreek",
         "Samford", "SilverPlains", "SturtPlains", "TiTreeEast", "Tumbarumba", "WallabyCreek", "Warra",
         "Whroo", "WombatStateForest", "Yanco"]
for site in sites:
    sp = os.path.join(rp, site, "Data", "Portal")
    op = os.path.join(rp, site, "Data", "Processed")
    if not os.path.isdir(sp):
        msg = sp + " , skipping site ..."
        logger.warning(msg)
        continue
    files = sorted([f for f in os.listdir(sp) if ("L3" in f and ".nc" in f)])
    #files = sorted([f for f in os.listdir(sp) if ("_L3.nc" in f)])
    if len(files) == 0:
        msg = "No files found in " + sp + " , skipping ..."
        logger.error(msg)
        continue
    for fn in files:
        ifp = os.path.join(sp, fn)
        msg = "Converting " + fn
        logger.info(msg)
        std["Files"]["in_filename"] = ifp
        # read the input file
        ds1 = pfp_io.NetCDFRead(ifp, update=False)
        # update the variable names
        change_variable_names(std, ds1)
        # make sure there are Ws and Wd series
        copy_ws_wd(ds1)
        # make sure we have all the variables we want ...
        ds2 = include_variables(std, ds1)
        # ... but not the ones we don't
        exclude_variables(std, ds2)
        # update the global attributes
        change_global_attributes(std, ds2)
        # update the variable attributes
        change_variable_attributes(std, ds2)
        # Fc single point storage
        consistent_Fco2_storage(std, ds2, site)
        ofp = os.path.join(op, fn)
        #nf = pfp_io.nc_open_write(ofp)
        #pfp_io.nc_write_series(nf, ds2)
        pfp_io.NetCDFWrite(ofp, ds2)
        # run cfchecker on the netCDF file
        cmd = ["cfchecks", "-v 1.8", ofp]
        #subprocess.run(cmd, stdout=cfchecker_file)
cfchecker_file.close()
# parse the cfchecker output file and write separate files for errors and warnings
error_file = open(error_file_name, "w")
warning_file = open(warning_file_name, "w")
with open(cfchecker_file_name) as f:
    for line in f:
        if "CHECKING NetCDF FILE" in line:
            error_file.write(line)
            warning_file.write(line)
            continue
        if "Checking variable:" in line:
            parts = line.split(" ")
            variable_name = parts[2].rstrip()
            continue
        if "ERROR:" in line:
            line_out = " " + variable_name + " " + line
            error_file.write(line_out)
            continue
        if "WARN:" in line:
            line_out = " " + variable_name + " " + line
            warning_file.write(line_out)
            continue
error_file.close()
warning_file.close()
