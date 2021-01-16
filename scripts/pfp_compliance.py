# standard modules
import copy
import datetime
import inspect
import logging
import ntpath
import os
import platform
import sys
import time
import traceback
# 3rd party modules
from configobj import ConfigObj
import numpy
from PyQt5 import QtWidgets
import scipy.stats
#import timezonefinder
import xlrd
# PFP modules
from scripts import constants as c
from scripts import pfp_cfg
from scripts import pfp_func_units
from scripts import pfp_func_stats
from scripts import pfp_gui
from scripts import pfp_io
from scripts import pfp_ts
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def change_variable_names(cfg, ds):
    """
    Purpose:
     Change variable names to the new (October 2018) scheme.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # get a list of potential mappings
    renames_exact = list(cfg["rename_exact"].keys())
    # loop over the variables in the data structure
    labels = list(ds.series.keys())
    for label in labels:
        if label in renames_exact:
            new_name = cfg["rename_exact"][label]
            ds.series[new_name] = ds.series.pop(label)
    # rename pattern matches
    renames_pattern = list(cfg["rename_pattern"].keys())
    labels = list(ds.series.keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = cfg["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label in labels:
            if ((label[:llen] == rp) and
                (label[:klen] != srp)):
                new_name = label.replace(label[:llen], srp)
                ds.series[new_name] = ds.series.pop(label)
    return

def change_variable_units(cfg, ds):
    """
    Purpose:
     Change variable units to the new (May 2020) scheme.
    Usage:
    Author: PRI
    Date: June 2020
    """
    labels = list(ds.series.keys())
    # coerce units into a standard form
    old_units = list(cfg["units_map"].keys())
    new_units = [cfg["units_map"][o] for o in old_units]
    ok_units = list(set(old_units + new_units))
    for label in labels:
        if label in ["DateTime", "time"]:
            continue
        units = ds.series[label]["Attr"]["units"]
        if len(units) == 0:
            continue
        if units not in ok_units:
            #msg = " Unrecognised units " + units + " for variable " + label
            #logger.warning(msg)
            continue
        if units in old_units:
            ds.series[label]["Attr"]["units"] = cfg["units_map"][units]
    return

def CheckExcelWorkbook(l1_info):
    """
    Purpose:
     Check that the requested sheets and variables are in the Excel workbook.
    Usage:
    Side effects:
    Author: PRI
    Date: February 2020
    """
    l1ire = l1_info["read_excel"]
    # get the labels of variables in the netCDF file
    # need to add check for duplicate variables on multiple sheets
    nc_labels = list(l1ire["Variables"].keys())
    l1ire["xl_sheets"] = {}
    # check the requested sheets are present and get a list of variables for each sheet
    msg = " Reading Excel workbook " + os.path.basename(l1ire["Files"]["in_filename"])
    logger.info(msg)
    xl_book = xlrd.open_workbook(l1ire["Files"]["file_name"], on_demand=True)
    # put the Excel workbook datemode into the global attributes
    l1ire["Global"]["xl_datemode"] = str(xl_book.datemode)
    xl_sheets_present = xl_book.sheet_names()
    for nc_label in nc_labels:
        if "xl" in list(l1ire["Variables"][nc_label].keys()):
            xl_sheet = l1ire["Variables"][nc_label]["xl"]["sheet"]
            xl_label = l1ire["Variables"][nc_label]["xl"]["name"]
            if xl_sheet not in xl_sheets_present:
                msg = " Sheet " + xl_sheet + " (" + xl_label + ") not found in workbook, skipping "
                msg += nc_label + " ..."
                logger.warning(msg)
                del l1ire["Variables"][nc_label]
                continue
            if xl_sheet not in list(l1ire["xl_sheets"].keys()):
                l1ire["xl_sheets"][xl_sheet] = {"DateTime": "", "xl_labels":{}}
            l1ire["xl_sheets"][xl_sheet]["xl_labels"][xl_label] = nc_label
    # check the requested variables are on the specified sheets
    xl_data = {}
    for xl_sheet in list(l1ire["xl_sheets"].keys()):
        xl_data[xl_sheet] = xl_book.sheet_by_name(xl_sheet)
        headers = xl_data[xl_sheet].row_values(l1ire["Files"]["in_headerrow"])
        for xl_label in list(l1ire["xl_sheets"][xl_sheet]["xl_labels"].keys()):
            if xl_label not in headers:
                msg = " Variable " + xl_label + " not found on sheet " + xl_sheet + ", skipping ..."
                logger.warning(msg)
                del l1ire["xl_sheets"][xl_sheet]["xl_labels"][xl_label]
    # now we know what variables on which sheets have been requested and are present
    # find the timestamp label for each sheet
    fdr = int(l1ire["Files"]["in_firstdatarow"])
    for xl_sheet in list(l1ire["xl_sheets"].keys()):
        got_timestamp = False
        ldr = int(xl_data[xl_sheet].nrows)
        xl_labels = xl_data[xl_sheet].row_values(l1ire["Files"]["in_headerrow"])
        for xl_label in xl_labels:
            col = xl_labels.index(xl_label)
            types = numpy.array(xl_data[xl_sheet].col_types(col)[fdr:ldr])
            mode = scipy.stats.mode(types)
            if mode[0][0] == xlrd.XL_CELL_DATE and 100*mode[1][0]//len(types) > 75:
                got_timestamp = True
                l1ire["xl_sheets"][xl_sheet]["DateTime"] = xl_label
                if xl_label in l1ire["xl_sheets"][xl_sheet]["xl_labels"]:
                    del l1ire["xl_sheets"][xl_sheet]["xl_labels"][xl_label]
                break
        if not got_timestamp:
            msg = " Time stamp not found for sheet " + xl_sheet +", skipping ..."
            logger.warning(msg)
            del l1ire["xl_sheets"][xl_sheet]

    return xl_data

def check_executables():
    # check for the executable files required
    if platform.system() == "Windows":
        executable_extension = ".exe"
    else:
        executable_extension = ""
    executable_names = ["solo/bin/sofm", "solo/bin/solo", "solo/bin/seqsolo",
                        "mds/bin/gf_mds", "mpt/bin/ustar_mp"]
    missing_executable = False
    for executable_name in executable_names:
        executable_name = executable_name + executable_extension
        if not os.path.isfile(executable_name):
            missing_executable = True
    if missing_executable:
        msg = "One or more of the required executable files could not be found.\n"
        msg = msg + "If you are running on Windows, clone the repository again.\n"
        msg = msg + "If you are running on OSX or Linux, use the make_nix script\n"
        msg = msg + "to compile the executables."
        result = pfp_gui.MsgBox_Quit(msg, title="Critical")
    return

def check_status_ok(ds, info):
    if info["status"]["value"] != 0:
        ds.returncodes["value"] = info["status"]["value"]
        ds.returncodes["message"] = info["status"]["message"]
        logger.error(info["status"]["message"])
        return False
    else:
        return True

def climatology_update_controlfile(cfg):
    """
    Purpose:
     Parse the climatology control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.climatology_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: September 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        strip_list = ['"', "'", "[", "]"]
        for key1 in cfg:
            if key1 in ["level"]:
                cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
            elif key1 in ["Files"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    cfg[key1][key2] = parse_cfg_values(key2, cfg2, strip_list)
            elif key1 in ["Variables"]:
                for key2 in cfg[key1]:
                    for key3 in cfg[key1][key2]:
                        cfg3 = cfg[key1][key2][key3]
                        cfg[key1][key2][key3] = parse_cfg_values(key3, cfg3, strip_list)
            else:
                del cfg[key1]
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    except Exception:
        ok = False
        msg = " An error occurred while updating the climatology control file syntax"
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def cpd_mchugh_update_controlfile(cfg):
    """
    Purpose:
     Parse the CPD (McHugh) control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.cpd_mchugh_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: September 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        strip_list = ['"', "'", "[", "]"]
        for key1 in cfg:
            if key1 in ["level"]:
                cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
            elif key1 in ["Files", "Options"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    cfg[key1][key2] = parse_cfg_values(key2, cfg2, strip_list)
            elif key1 in ["Variables"]:
                for key2 in cfg[key1]:
                    for key3 in cfg[key1][key2]:
                        cfg3 = cfg[key1][key2][key3]
                        cfg[key1][key2][key3] = parse_cfg_values(key3, cfg3, strip_list)
            else:
                del cfg[key1]
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    except Exception:
        ok = False
        msg = " An error occurred while updating the CPD (McHugh) control file syntax"
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def cpd_barr_update_controlfile(cfg):
    """
    Purpose:
     Parse the CPD (Barr) control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.cpd_barr_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: September 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        strip_list = ['"', "'", "[", "]"]
        for key1 in cfg:
            if key1 in ["level"]:
                cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
            elif key1 in ["Files", "Options"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    cfg[key1][key2] = parse_cfg_values(key2, cfg2, strip_list)
            elif key1 in ["Variables"]:
                for key2 in cfg[key1]:
                    for key3 in cfg[key1][key2]:
                        cfg3 = cfg[key1][key2][key3]
                        cfg[key1][key2][key3] = parse_cfg_values(key3, cfg3, strip_list)
            else:
                del cfg[key1]
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    except Exception:
        ok = False
        msg = " An error occurred while updating the CPD (Barr) control file syntax"
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def mpt_update_controlfile(cfg):
    """
    Purpose:
     Parse the MPT control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.mpt_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: September 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        strip_list = ['"', "'", "[", "]"]
        for key1 in cfg:
            if key1 in ["level"]:
                cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
            elif key1 in ["Files", "Options"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    cfg[key1][key2] = parse_cfg_values(key2, cfg2, strip_list)
            elif key1 in ["Variables"]:
                for key2 in cfg[key1]:
                    for key3 in cfg[key1][key2]:
                        cfg3 = cfg[key1][key2][key3]
                        cfg[key1][key2][key3] = parse_cfg_values(key3, cfg3, strip_list)
            else:
                del cfg[key1]
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    except Exception:
        ok = False
        msg = " An error occurred while updating the MPT control file syntax"
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def concatenate_update_controlfile(cfg):
    """
    Purpose:
     Parse the concatenate control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.concatenate_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: February 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        strip_list = ['"', "'", "[", "]"]
        for key1 in cfg:
            if key1 in ["level"]:
                cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
            elif key1 in ["Options"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    if key2 in ["DoFingerprints", "FixTimeStepMethod", "MaxGapInterpolate",
                                "NumberOfDimensions", "SeriesToCheck", "Truncate",
                                "TruncateThreshold"]:
                        cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                    else:
                        del cfg[key1][key2]
            elif key1 in ["Files"]:
                for key2 in cfg[key1]:
                    cfg2 = cfg[key1][key2]
                    for key3 in cfg2:
                        cfg3 = cfg[key1][key2][key3]
                        cfg3 = parse_cfg_values(key3, cfg3, strip_list)
            else:
                del cfg[key1]
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    except Exception:
        ok = False
        msg = " An error occurred while updating the concatenate control file syntax"
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def consistent_Fco2_storage(ds, file_name):
    """
    Purpose:
     Make the various incarnations of single point Fco2 storage consistent.
    Author: PRI
    Date: November 2019
    """
    labels = list(ds.series.keys())
    if "CO2" not in labels:
        # do nothing if CO2 concentration does not exist
        pass
    elif "Fco2_single" in labels:
        # do nothing if Fco2_single already exists
        pass
    elif "Fco2_storage" in labels:
        # Fc_single may be called Fc_storage in earlier files
        level = ds.globalattributes["nc_level"]
        descr = "description_" + level
        variable = pfp_utils.GetVariable(ds, "Fco2_storage")
        if "using single point CO2 measurement" in variable["Attr"][descr]:
            variable["Label"] = "Fco2_single"
            pfp_utils.CreateVariable(ds, variable)
            pfp_utils.DeleteVariable(ds, "Fco2_storage")
    else:
        # neither Fco2_single nor Fco2_storage exist, try to calculate
        # check to see if the measurement height is defined
        CO2 = pfp_utils.GetVariable(ds, "CO2")
        zms = pfp_utils.get_number_from_heightstring(CO2["Attr"]["height"])
        if zms == None:
            while zms == None:
                text, ok = QtWidgets.QInputDialog.getText(None, file_name,
                                                          "Enter CO2 measuement height in metres",
                                                          QtWidgets.QLineEdit.Normal,"")
                zms = pfp_utils.get_number_from_heightstring(text)
        # update the CO2 variable attribute
        CO2["Attr"]["height"] = zms
        pfp_utils.CreateVariable(ds, CO2)
        # calculate single point Fc storage term
        cf = {"Options": {"zms": zms}}
        pfp_ts.CalculateFco2StorageSinglePoint(cf, ds)
        # convert Fco2_single from mg/m^2/s to umol/m^2/s
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
    # get a list of the series
    series_list = sorted(list(ds.series.keys()))
    if "Wd" not in series_list:
        if "Wd_SONIC_Av" in series_list:
            ds.series["Wd"] = copy.deepcopy(ds.series["Wd_SONIC_Av"])
            ds.series["Wd"]["Attr"]["long_name"] = "Wind direction (copied from Wd_SONIC_Av)"
    if "Ws" not in series_list:
        if "Ws_SONIC_Av" in series_list:
            ds.series["Ws"] = copy.deepcopy(ds.series["Ws_SONIC_Av"])
            ds.series["Ws"]["Attr"]["long_name"] = "Wind speed (copied from Ws_SONIC_Av)"
    return

def ParseConcatenateControlFile(cf):
    """
    Purpose:
     Make the concatenate information dictionary
    Usage:
    Side effects:
    Author: PRI
    Date: August 2019
    """
    info = {}
    info["NetCDFConcatenate"] = {"OK": True}
    inc = info["NetCDFConcatenate"]
    # check the control file has a Files section
    if "Files" not in cf:
        msg = " Files section missing from control file"
        logger.error(msg)
        inc["OK"] = False
        return info
    # check the [Files] section contains an [Out] section and an [In] section
    for item in ["Out", "In"]:
        if item not in cf["Files"]:
            msg = " " + item + " subsection missing from Files section"
            logger.error(msg)
            inc["OK"] = False
            return info
    # check the [In] section contains at least 1 entry
    if len(list(cf["Files"]["In"].keys())) < 2:
        msg = " Less than 2 input files specified"
        logger.error(msg)
        inc["OK"] = False
        return info
    # get a list of the input file names
    inc["in_file_names"] = []
    for key in sorted(list(cf["Files"]["In"].keys())):
        file_name = cf["Files"]["In"][key]
        if os.path.isfile(file_name):
            inc["in_file_names"].append(file_name)
        else:
            msg = " File not found (" + ntpath.basename(file_name) + ")"
            logger.warning(msg)
    # get the output file name
    if "ncFileName" not in cf["Files"]["Out"]:
        msg = " No ncFileName key in Out subsection of Files section"
        logger.error(msg)
        inc["OK"] = False
        return info
    inc["out_file_name"] = cf["Files"]["Out"]["ncFileName"]
    # check the output path exists, create if it doesn't
    file_path, file_name = os.path.split(inc["out_file_name"])
    if not os.path.isdir(file_path):
        os.makedirs(file_path)
    # work through the choices in the [Options] section
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "NumberOfDimensions", default=3)
    inc["NumberOfDimensions"] = int(opt)
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "MaxGapInterpolate", default=0)
    inc["MaxGapInterpolate"] = int(opt)
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "FixTimeStepMethod", default="round")
    inc["FixTimeStepMethod"] = str(opt)
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "Truncate", default="No")
    inc["Truncate"] = str(opt)
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "TruncateThreshold", default=50)
    inc["TruncateThreshold"] = float(opt)
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "SeriesToCheck", default="all")
    inc["SeriesToCheck"] = pfp_utils.csv_string_to_list(opt)
    # now add the bits and pieces
    inc["start_date"] = []
    inc["end_date"] = []
    inc["chrono_files"] = []
    inc["labels"] = []
    inc["attributes"] = ["group_name", "height", "instrument", "long_name",
                         "standard_name", "units", "valid_range"]
    # add key for updating netCDF files
    stdname = os.path.join("controlfiles", "standard", "nc_cleanup.txt")
    info["NetCDFUpdate"] = pfp_io.get_controlfilecontents(stdname)
    # add key for suppressing output of intermediate variables e.g. Cpd etc
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "KeepIntermediateSeries", default="No")
    info["RemoveIntermediateSeries"] = {"KeepIntermediateSeries": opt, "not_output": []}
    return info

def ParseL1ControlFile(cf):
    """
    Purpose:
     Check the contents of the Li control file.
     If the L1 control file contents are OK, return with the required information.
     If the L1 control file contents are not OK, return with an error message.
    """
    logger.info(" Parsing the L1 control file")
    # create the settings dictionary
    l1_info = {"status": {"value": 0, "message": "OK"},
              "read_excel": {}}
    l1ire = l1_info["read_excel"]
    # copy the files section from the control file
    l1ire["Files"] = copy.deepcopy(cf["Files"])
    l1ire["Files"]["file_name"] = os.path.join(cf["Files"]["file_path"], cf["Files"]["in_filename"])
    l1ire["Files"]["in_headerrow"] = int(cf["Files"]["in_headerrow"])- 1
    l1ire["Files"]["in_firstdatarow"] = int(cf["Files"]["in_firstdatarow"]) - 1
    # get the global attributes
    l1ire["Global"] = copy.deepcopy(cf["Global"])
    # get the variables
    l1ire["Variables"] = copy.deepcopy(cf["Variables"])
    # checks of the 'Variables' sections would go here
    for label in list(l1ire["Variables"].keys()):
        # check the 'xl' subsection is present
        ok = True
        if "xl" not in list(l1ire["Variables"][label].keys()):
            # the xl subsection can be missing if Function is present
            if "Function" in list(l1ire["Variables"][label].keys()):
                pass
            else:
                msg = " Skipping " + label + " (subsections 'xl' or 'Function' not found)"
                logger.warning(msg)
                ok = False
        # check the 'Attr' subsection is present
        if "Attr" not in list(l1ire["Variables"][label].keys()):
            msg = " Skipping " + label + " (subsection 'Attr' not found)"
            logger.warning(msg)
            ok = False
        if not ok:
            del l1ire["Variables"][label]
            continue
        # check 'sheet' and 'name' are in the 'xl' subsection
        ok = True
        if "xl" in list(l1ire["Variables"][label].keys()):
            for item in ["sheet", "name"]:
                if item not in list(l1ire["Variables"][label]["xl"].keys()):
                    msg = " Skipping " + label + " ('" + item + "' not found in 'xl' subsection)"
                    logger.warning(msg)
                    ok = False
        # check the 'Function' subsection
        elif "Function" in list(l1ire["Variables"][label].keys()):
            # check 'func' is in the 'Function' subsection
            if "func" not in list(l1ire["Variables"][label]["Function"].keys()):
                msg = " Skipping " + label + " ('func' not found in 'Function' subsection)"
                logger.warning(msg)
                ok = False
            # check the function name and arguments
            else:
                # get a list of function names in pfp_func_units
                implemented_func_units = [name for name,data in inspect.getmembers(pfp_func_units,inspect.isfunction)]
                implemented_func_stats = [name for name,data in inspect.getmembers(pfp_func_stats,inspect.isfunction)]
                implemented_functions = implemented_func_units + implemented_func_stats
                function_string = l1ire["Variables"][label]["Function"]["func"]
                function_name = function_string[:function_string.index("(")]
                # check the function name is implemented
                if function_name not in implemented_functions:
                    msg = " Skipping " + label + " (function " + function_name + " not implemented)"
                    logger.warning(msg)
                    ok = False
                # check the arguments are being read in
                else:
                    function_args = function_string[function_string.index("(")+1:-1].split(",")
                    nargs = len(function_args)
                    if function_name in ["Linear"]:
                        nargs = 1
                    for item in function_args[:nargs]:
                        if item not in list(l1ire["Variables"].keys()):
                            msg = " Skipping " + label + " (function argument '" + item + "' not found)"
                            logger.warning(msg)
                            ok = False
                        else:
                            pass
        # we should never get here
        else:
            msg = " These are not the droids you are looking for!"
            logger.error(msg)
            ok = False
        if not ok:
            del l1ire["Variables"][label]
            continue
        # check the 'Attr' subsection
        ok = True
        for item in ["long_name", "units"]:
            if item not in list(l1ire["Variables"][label]["Attr"].keys()):
                msg = " Skipping " + label + " (subsection '" + item + "' not found)"
                logger.warning(msg)
                ok = False
        if not ok:
            del l1ire["Variables"][label]
            continue
    return l1_info

def ParseL3ControlFile(cf, ds):
    """
    Purpose:
     Parse the L3 control file and return contents in the l3_info dictionary.
    Usage:
    Side effects:
    Author: PRI
    Date: August 2019
    """
    ds.returncodes["message"] = "OK"
    ds.returncodes["value"] = 0
    l3_info = {"CO2": {}, "Fco2": {}, "status": {"value": 0, "message": "OK"}}
    # add key for suppressing output of intermediate variables e.g. Cpd etc
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "KeepIntermediateSeries", default="No")
    l3_info["RemoveIntermediateSeries"] = {"KeepIntermediateSeries": opt, "not_output": []}
    # find out what label is used for CO2
    if "CO2" in list(cf["Variables"].keys()):
        l3_info["CO2"]["label"] = "CO2"
    else:
        msg = " Label for CO2 not found in control file"
        l3_info["status"]["value"] = 1
        l3_info["status"]["message"] = msg
        return l3_info
    # get the height of the CO2 measurement
    got_zms = False
    labels = list(ds.series.keys())
    CO2_label = l3_info["CO2"]["label"]
    # try and get the height from the CO2 variable
    if CO2_label in labels:
        # get height from attributes if the CO2 variable is already in the data structure
        try:
            CO2 = pfp_utils.GetVariable(ds, CO2_label)
            zms = float(pfp_utils.strip_non_numeric(CO2["Attr"]["height"]))
            got_zms = True
        except:
            pass
    if not got_zms and "MergeSeries" in list(cf["Variables"][CO2_label].keys()):
        # get the height from the variables listed in MergeSeries
        try:
            source = cf["Variables"][CO2_label]["MergeSeries"]["source"]
            source = pfp_utils.convert_csv_string_to_list(source)
            for item in source:
                var = pfp_utils.GetVariable(ds, item)
                zms = float(pfp_utils.strip_non_numeric(var["Attr"]["height"]))
                got_zms = True
                break
        except:
            pass
    if not got_zms and "AverageSeries" in list(cf["Variables"][CO2_label].keys()):
        # get the height from the variables listed in AverageSeries
        try:
            source = cf["Variables"][CO2_label]["AverageSeries"]["source"]
            source = pfp_utils.convert_csv_string_to_list(source)
            for item in source:
                var = pfp_utils.GetVariable(ds, item)
                zms = float(pfp_utils.strip_non_numeric(var["Attr"]["height"]))
                got_zms = True
                break
        except:
            pass
    if not got_zms and "tower_height" in list(ds.globalattributes.keys()):
        try:
            zms = float(pfp_utils.strip_non_numeric(ds.globalattributes["tower_height"]))
            got_zms = True
        except:
            pass
    if not got_zms and pfp_utils.cfkeycheck(cf, Base="Options", ThisOne="zms"):
        try:
            zms = float(pfp_utils.strip_non_numeric(cf["Options"]["zms"]))
            got_zms = True
        except:
            pass
    if got_zms:
        l3_info["CO2"]["height"] = zms
    else:
        msg = " Unable to find height for CO2 (" + CO2_label + ") measurement"
        l3_info["status"]["value"] = 1
        l3_info["status"]["message"] = msg
        return l3_info
    # get a list of Fco2 variables to be merged
    cfv = cf["Variables"]
    merge_list = [l for l in list(cfv.keys()) if l[0:4] == "Fco2" and "MergeSeries" in list(cfv[l].keys())]
    average_list = [l for l in list(cfv.keys()) if l[0:4] == "Fco2" and "AverageSeries" in list(cfv[l].keys())]
    l3_info["Fco2"]["combine_list"] = merge_list + average_list
    return l3_info

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
        if not isinstance(value, basestring):
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

def remove_variables(cfg, ds):
    """
    Purpose:
     Remove deprecated variables from a netCDF file.
    Usage:
    Author: PRI
    Date: October 2018
    """
    remove_list = [v for v in list(cfg["Variables"].keys()) if "remove" in cfg["Variables"][v]]
    series_list = sorted(list(ds.series.keys()))
    for label in series_list:
        if label in remove_list:
            ds.series.pop(label)
    return

def change_global_attributes(cfg, ds):
    """
    Purpose:
     Clean up the global attributes.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # check site_name is in ds.globalattributes
    gattr_list = list(ds.globalattributes.keys())
    if "site_name" not in gattr_list:
        logger.warning("Global attributes: site_name not found")
    # check latitude and longitude are in ds.globalattributes
    if "latitude" not in gattr_list:
        logger.warning("Global attributes: latitude not found")
    else:
        lat_string = str(ds.globalattributes["latitude"])
        if len(lat_string) == 0:
            logger.warning("Global attributes: latitude empty")
        else:
            lat = pfp_utils.convert_anglestring(lat_string)
        ds.globalattributes["latitude"] = str(lat)
    if "longitude" not in gattr_list:
        logger.warning("Global attributes: longitude not found")
    else:
        lon_string = str(ds.globalattributes["longitude"])
        if len(lon_string) == 0:
            logger.warning("Global attributes: longitude empty")
        else:
            lon = pfp_utils.convert_anglestring(lon_string)
        ds.globalattributes["longitude"] = str(lon)
    # check to see if there there is a time_zone global attribute
    gattr_list = list(ds.globalattributes.keys())
    #if not "time_zone" in gattr_list:
        ## get the site name
        #site_name = ds.globalattributes["site_name"]
        #sn = site_name.replace(" ","").replace(",","").lower()
        ## first, see if the site is in constants.tz_dict
        #if sn in list(c.tz_dict.keys()):
            #ds.globalattributes["time_zone"] = c.tz_dict[sn]
        #else:
            #if "latitude" in gattr_list and "longitude" in gattr_list:
                #lat = float(ds.globalattributes["latitude"])
                #lon = float(ds.globalattributes["longitude"])
                #if lat != -9999 and lon != -9999:
                    #tf = timezonefinder.TimezoneFinder()
                    #tz = tf.timezone_at(lng=lon, lat=lat)
                    #ds.globalattributes["time_zone"] = tz
                #else:
                    #logger.warning("Global attributes: unable to define time zone")
                    #ds.globalattributes["time_zone"] = ""
    # add or change global attributes as required
    gattr_list = sorted(list(cfg["Global"].keys()))
    for gattr in gattr_list:
        ds.globalattributes[gattr] = cfg["Global"][gattr]
    # remove deprecated global attributes
    flag_list = [g for g in list(ds.globalattributes.keys()) if "Flag" in g]
    others_list = ["end_datetime", "start_datetime", "Functions", "doi"]
    remove_list = others_list + flag_list
    for gattr in list(ds.globalattributes.keys()):
        if gattr in remove_list:
            ds.globalattributes.pop(gattr)
    # rename global attributes
    rename_dict = {"EPDversion":"PythonVersion", "elevation":"altitude"}
    for item in rename_dict:
        if item in list(ds.globalattributes.keys()):
            new_key = rename_dict[item]
            ds.globalattributes[new_key] = ds.globalattributes.pop(item)
    return

def change_variable_attributes(cfg, ds):
    """
    Purpose:
     Clean up the variable attributes.
    Usage:
    Author: PRI
    Date: November 2018
    """
    # rename existing long_name to description, introduce a
    # consistent long_name attribute and introduce the group_name
    # attribute
    vattr_list = list(cfg["variable_attributes"].keys())
    series_list = list(ds.series.keys())
    if "nc_level" not in ds.globalattributes:
        ds.globalattributes["nc_level"] = "L1"
    descr = "description_" + ds.globalattributes["nc_level"]
    for label in series_list:
        variable = pfp_utils.GetVariable(ds, label)
        variable["Attr"][descr] = copy.deepcopy(variable["Attr"]["long_name"])
        for item in vattr_list:
            if label[:len(item)] == item:
                for key in list(cfg["variable_attributes"][item].keys()):
                    if key != "units":
                        variable["Attr"][key] = cfg["variable_attributes"][item][key]
        pfp_utils.CreateVariable(ds, variable)
    # parse variable attributes to new format, remove deprecated variable attributes
    # and fix valid_range == "-1e+35,1e+35"
    tmp = cfg["deprecated"]["attributes"]
    deprecated_attributes = pfp_cfg.cfg_string_to_list(tmp)
    tmp = cfg["deprecated"]["values"]
    deprecated_values = pfp_cfg.cfg_string_to_list(tmp)
    series_list = list(ds.series.keys())
    for label in series_list:
        variable = pfp_utils.GetVariable(ds, label)
        # parse variable attributes to new format
        variable["Attr"] = parse_variable_attributes(variable["Attr"])
        # remove deprecated variable attributes
        for vattr in deprecated_attributes:
            if vattr in list(variable["Attr"].keys()):
                del variable["Attr"][vattr]
        # remove attributes with deprecated values
        for vattr in list(variable["Attr"].keys()):
            if variable["Attr"][vattr] in deprecated_values:
                del variable["Attr"][vattr]
        # remove attributes with empty values
        vattrs_essential = pfp_cfg.cfg_string_to_list(cfg["essential"]["attributes"])
        for vattr in list(variable["Attr"].keys()):
            if ((len(str(variable["Attr"][vattr])) == 0) and (vattr not in vattrs_essential)):
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
    return

def error_handler(info, msg, val):
    """
    Purpose:
     Handle errors from ParseL1ControlFile().
    """
    logger.error(msg)
    info["status"]["message"] = msg
    info["status"]["value"] = val
    return

def exclude_variables(cfg, ds):
    """
    Purpose:
     Remove deprecated variables from a netCDF file.
    Usage:
    Author: PRI
    Date: October 2018
    """
    series_list = sorted(list(ds.series.keys()))
    var_list = [v for v in list(cfg["exclude"].keys())]
    flag_list = [v+"_QCFlag" for v in var_list if v+"_QCFlag" in series_list]
    remove_list = var_list + flag_list
    for label in series_list:
        if label in remove_list:
            ds.series.pop(label)
    return

def include_variables(cfg, ds_in):
    """
    Purpose:
     Only pick variables that match the specified string for the length
     of the specified string.
    Usage:
    Author: PRI
    Date: November 2018
    """
    # get a new data structure
    ds_out = pfp_io.DataStructure()
    # copy the global attributes
    for gattr in ds_in.globalattributes:
        ds_out.globalattributes[gattr] = ds_in.globalattributes[gattr]
    # loop over variables to be included
    include_list = list(cfg["include"].keys())
    series_list = list(ds_in.series.keys())
    for item in include_list:
        for label in series_list:
            if label[0:len(item)] == item:
                ds_out.series[label] = ds_in.series[label]
    return ds_out

def check_batch_controlfile(cfg):
    """
    Purpose:
     Check the L1 control file to make sure it contains all information
     needed to run L1 and that all information is correct.
    Usage:
    Side effects:
    Author: PRI
    Date: June 2020
    """
    ok = True
    return ok

def check_l1_controlfile(cfg):
    """
    Purpose:
     Check the L1 control file to make sure it contains all information
     needed to run L1 and that all information is correct.
    Usage:
    Side effects:
    Author: PRI
    Date: June 2020
    """
    ok = True
    # check the expected sections are in the control file
    for item in ["Files", "Global", "Variables"]:
        if item not in list(cfg.keys()):
            msg = " Section '" + item + "' not in control file"
            logger.error(msg)
            ok = False
    # check the directory exists
    file_path = pfp_utils.get_keyvaluefromcf(cfg, ["Files"], "file_path", default=".")
    if not os.path.isdir(file_path):
        msg = " Input directory " + file_path + " not found"
        error_message(msg, mode="correct")
        ok = False
    # get the input file and check it exists
    xl_file_name = pfp_io.get_infilenamefromcf(cfg)
    if not os.path.isfile(xl_file_name):
        msg = " Input file " + os.path.basename(xl_file_name) + " not found"
        error_message(msg, mode="correct")
        ok = False
    # check the header row entry are numbers
    try:
        opt = pfp_utils.get_keyvaluefromcf(cfg, ["Files"], "in_headerrow", default=2)
        opt = int(opt) - 1
    except ValueError:
        msg = " In the Files section of the control file, in_headerrow is not a number"
        error_message(msg, mode="correct")
        ok = False
    # check the first data row entry are numbers
    try:
        opt = pfp_utils.get_keyvaluefromcf(cfg, ["Files"], "in_firstdatarow", default=5)
        opt = int(opt) - 1
    except ValueError:
        msg = " In the Files section of the control file, in_firstdatarow is not a number"
        error_message(msg, mode="correct")
        ok = False
    # *************************
    # *** global attributes ***
    # *************************
    # check time step is present and makes sense
    try:
        ts = int(cfg["Global"]["time_step"])
    except ValueError:
        msg = " Global attribute 'time_step' is not a number"
        error_message(msg, mode="correct")
        ok = False
    if ts not in [15, 20, 30, 60]:
        msg = " Global attribute 'time_step' must be 15, 20, 30 or 60"
        error_message(msg, mode="correct")
        ok = False
    # check latitude and longitude are present and make sense
    try:
        lat = float(cfg["Global"]["latitude"])
    except ValueError:
        msg = " Global attribute 'latitude' is not a number"
        error_message(msg, mode="correct")
        ok = False
    if lat < -90.0 or lat > 90.0:
        msg = "Global attribute 'latitude' must be between -90 and 90"
        error_message(msg, mode="correct")
        ok = False
    try:
        lon = float(cfg["Global"]["longitude"])
    except ValueError:
        msg = " Global attribute 'longitude' is not a number"
        error_message(msg, mode="correct")
        ok = False
    if lon < -180.0 or lat > 180.0:
        msg = " Global attribute 'longitude' must be between -180 and 180"
        error_message(msg, mode="correct")
        ok = False
    ## check the time zone
    #try:
        #tf = timezonefinder.TimezoneFinder()
        #tz_from_lat_lon = tf.timezone_at(lng=lon, lat=lat)
        #opt = pfp_utils.get_keyvaluefromcf(cfg, ["Global"], "time_zone", default=tz_from_lat_lon)
        #if opt.lower() != tz_from_lat_lon.lower():
            #msg = " Global attribute 'time_zone' inconsistent with latitude and longitude"
            #logger.warning(msg)
            #msg = "  " + opt + " replaced with " + tz_from_lat_lon
            #logger.warning(msg)
            #cfg["Global"]["time_zone"] = tz_from_lat_lon
            #file_name = os.path.basename(cfg.filename)
            #msg = " Updated 'time_zone' global attribute and saved control file " + file_name
            #logger.info(msg)
            #cfg.write()
    #except:
        #msg = " Error checking global attribute 'time_zone'"
        #error_message(msg)
        #ok = False
    # *******************
    # *** check units ***
    # *******************
    # the lists of units below need to be abstracted from the code and reconciled with
    # the contents of controlfiles/standard/cfg_update.txt
    units = {"co2": ["mg/m^3", "mmol/m^3", "umol/mol", "mg^2/m^6", "mmol^2/m^6",
                     "mg/m^2/s", "umol/m^2/s", "umol^2/mol^2"],
             "h2o": ["g/m^3", "kg/m^3", "mmol/m^3", "mmol/mol", "percent", "fraction",
                     "kg/kg", "g^2/m^6", "mmol^2/m^6", "mmol^2/mol^2"],
             "temperature": ["degC", "K", "degC^2", "K^2"],
             "pressure": ["Pa", "hPa", "kPa"],
             "soil": ["m^3/m^3", "dS/m"],
             "radiation": ["W/m^2", "umol/m^2/s", "mmol/m^2/s"],
             "covariance": ["g/m^2/s", "mg/m^2/s", "m.degC/s", "m.K/s", "m^2/s^2"],
             "flux": ["kg/m/s^2"],
             "precipitation": ["m", "mm"],
             "wind": ["m/s", "degrees", "m^2/s^2"],
             "heat": ["J/kg", "J/kg/K", "J/m^3/K"],
             "misc": ["V", "none", "ohms"]}
    ok_units = []
    for key in list(units.keys()):
        ok_units += units[key]
    ok_units = list(set(ok_units))
    for label in list(cfg["Variables"].keys()):
        var_units = cfg["Variables"][label]["Attr"]["units"]
        if var_units not in ok_units:
            msg = " Unrecognised units (" + var_units +") for variable " + label
            logger.error(msg)
            ok = False
    return ok

def error_message(msg, mode="correct"):
    logger.error("!!!")
    logger.error(msg)
    if mode == "correct":
        msg = " You can go back to the control file tab and correct the entry"
        logger.error(msg)
    logger.error("!!!")
    return

def l1_update_controlfile(cfg):
    """
    Purpose:
     Parse the L1 control file to update the syntax from earlier OFQC/PFP versions
     to the syntax used by this version.
    Usage:
     result = pfp_compliance.l1_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L1 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: February 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        cfg = l1_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred updating the L1 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        cfg = l1_update_cfg_variable_deprecate(cfg, std)
        cfg = l1_update_cfg_variable_names(cfg, std)
        cfg = l1_update_cfg_global_attributes(cfg, std)
        cfg = l1_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L1 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l1_update_cfg_global_attributes(cfg, std):
    """
    Purpose:
     Update the global attributes according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    # check for the essential global attributes
    essentials = ["latitude", "longitude", "site_name", "time_step", "time_zone"]
    for gattr in essentials:
        if gattr not in cfg["Global"]:
            cfg["Global"][gattr] = ""
    # remove deprecated global attributes
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["global"])
    for gattr in deprecated:
        if gattr in cfg["Global"]:
            cfg["Global"].pop(gattr)
    # rename global attributes
    renames = {"EPDversion":"PythonVersion", "elevation":"altitude"}
    for rename in list(renames.keys()):
        if rename in cfg["Global"]:
            cfg["Global"][renames[rename]] = cfg["Global"].pop(rename)
    # add or change global attributes as required
    gattrs = sorted(list(std["Global"].keys()))
    for gattr in gattrs:
        cfg["Global"][gattr] = std["Global"][gattr]
    return cfg

def l1_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L1 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: 9th May 2020, the day after my 64th birthday!
    """
    strip_list = ['"', "'", "[", "]"]
    if "level" not in list(cfg.keys()):
        cfg["level"] = "L1"
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files", "Global"]:
            for key2 in cfg[key1]:
                cfg[key1][key2] = parse_cfg_values(key2, cfg[key1][key2], strip_list)
        elif key1 in ["Variables"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    if key3 in ["xl", "csv", "Attr"]:
                        cfg3 = cfg[key1][key2][key3]
                        for key4 in cfg3:
                            # for keywords to lower case
                            if key4.lower() != key4:
                                cfg3[key4.lower()] = cfg3.pop(key4)
                            cfg[key1][key2][key3][key4.lower()] = parse_cfg_variables_value(key3, cfg3[key4.lower()])
        else:
            del cfg[key1]
    return cfg

def l1_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes.
    Usage:
    Comments: There must be a better way to do this that avoids the double
              "for" loop.
    Author: PRI
    Date: May 2020
    """
    # list of essential variable attributes
    vattrs_essential = pfp_cfg.cfg_string_to_list(std["essential"]["attributes"])
    vattrs_deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["attributes"])
    # list of standard attribute values
    labels_std = list(std["variable_attributes"].keys())
    labels_cfg = list(cfg["Variables"].keys())
    # add any essential variable attributes that are missing, deprecate those no longer used
    for label in labels_cfg:
        vattrs_cfg = list(cfg["Variables"][label]["Attr"].keys())
        for vattr in vattrs_essential:
            if vattr not in vattrs_cfg:
                cfg["Variables"][label]["Attr"][vattr] = "none"
        for vattr in vattrs_deprecated:
            if vattr in vattrs_cfg:
                del cfg["Variables"][label]["Attr"][vattr]
    # coerce units into a standard form
    old_units = list(std["units_map"].keys())
    new_units = [std["units_map"][o] for o in old_units]
    ok_units = list(set(old_units + new_units))
    for label_cfg in labels_cfg:
        cfg_units = cfg["Variables"][label_cfg]["Attr"]["units"]
        if cfg_units.lower() == "none":
            continue
        elif len(cfg_units) == 0:
            cfg["Variables"][label_cfg]["Attr"]["units"] = "none"
        elif cfg_units in old_units:
            cfg["Variables"][label_cfg]["Attr"]["units"] = std["units_map"][cfg_units]
        elif cfg_units not in ok_units:
            msg = " Unrecognised units " + cfg_units + " for variable " + label_cfg
            logger.warning(msg)
            continue
    # force some variable attributes to particular values
    for label_std in labels_std:
        # length of the label stub in the standard control file
        llen = len(label_std)
        # pointer to attributes in standard control file
        attr_std = std["variable_attributes"][label_std]
        # list of permitted units for variables that match this stub
        units_std = pfp_cfg.cfg_string_to_list(attr_std["units"])
        # loop over variables in the L1 control file
        labels_cfg = [l for l in list(cfg["Variables"].keys()) if l[:len(label_std)] == label_std]
        for label_cfg in labels_cfg:
            # pointer to attributes in user control file
            attr_cfg = cfg["Variables"][label_cfg]["Attr"]
            # units string given in the L1 control file
            units_cfg = attr_cfg["units"]
            if ((label_cfg[:llen] == label_std) and (units_cfg in units_std)):
                # the first letters and units match so update the long_name
                attr_cfg["long_name"] = attr_std["long_name"]
                attr_cfg["group_name"] = attr_std["group_name"]
                if "statistic_type" in list(attr_std.keys()):
                    attr_cfg["statistic_type"] = attr_std["statistic_type"]
                else:
                    attr_cfg["statistic_type"] = "average"
            elif (label_cfg[:llen] == label_std) and (label_cfg[-3:] == "_Sd"):
                attr_cfg["long_name"] = attr_std["long_name"]
                attr_cfg["group_name"] = attr_std["group_name"]
                if "statistic_type" in list(attr_std.keys()):
                    attr_cfg["statistic_type"] = attr_std["statistic_type"]
                else:
                    attr_cfg["statistic_type"] = "standard deviation"
            elif (label_cfg[:llen] == label_std) and (label_cfg[-3:] == "_Vr"):
                attr_cfg["long_name"] = attr_std["long_name"]
                attr_cfg["group_name"] = attr_std["group_name"]
                if "statistic_type" in list(attr_std.keys()):
                    attr_cfg["statistic_type"] = attr_std["statistic_type"]
                else:
                    attr_cfg["statistic_type"] = "variance"
            else:
                continue
    return cfg

def l1_update_cfg_variable_deprecate(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L1 control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["variables"])
    for label in deprecated:
        if label in cfg["Variables"]:
            cfg["Variables"].pop(label)
    return cfg

def l1_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    # dictionary of renamed variables
    renamed = {}
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # loop over the variables in the control file
    for label_cfg in list(cfg["Variables"].keys()):
        if label_cfg in renames_exact:
            new_name = std["rename_exact"][label_cfg]
            renamed[label_cfg] = new_name
            cfg["Variables"].rename(label_cfg, new_name)
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg["Variables"].keys()):
            if ((label_cfg[:llen] == rp) and (label_cfg[:klen] != srp)):
                if len(label_cfg) > len(rp):
                    if label_cfg[len(rp)] == "_":
                        new_name = label_cfg.replace(label_cfg[:llen], srp)
                        renamed[label_cfg] = new_name
                        cfg["Variables"].rename(label_cfg, new_name)
                    else:
                        # different variable name, leave it alone
                        pass
                else:
                    new_name = label_cfg.replace(label_cfg[:llen], srp)
                    renamed[label_cfg] = new_name
                    cfg["Variables"].rename(label_cfg, new_name)
    # do any functions
    for label_cfg in list(cfg["Variables"].keys()):
        if "Function" in cfg["Variables"][label_cfg]:
            func_str = cfg["Variables"][label_cfg]["Function"]["func"]
            for old_label in list(renamed.keys()):
                if old_label in func_str:
                    func_str = func_str.replace(old_label, renamed[old_label])
                    cfg["Variables"][label_cfg]["Function"]["func"] = func_str
    return cfg

def l2_update_controlfile(cfg):
    """
    Purpose:
     Parse the L2 control file to update the syntax from earlier OFQC/PFP versions
     to the syntax used by this version.
    Usage:
     result = pfp_compliance.l2_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L2 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: February 2020
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        cfg = l2_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L2 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        cfg = l2_update_cfg_variable_deprecate(cfg, std)
        cfg = l2_update_cfg_variable_names(cfg, std)
        cfg = l2_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L2 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l2_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L2 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: May 2020
    """
    strip_list = ['"', "'", "[", "]"]
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            for key2 in cfg[key1]:
                if key2 in ["irga_type", "SONIC_Check", "IRGA_Check"]:
                    cfg2 = cfg[key1][key2]
                    cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                    cfg[key1][key2] = cfg2
                else:
                    del cfg[key1][key2]
        elif key1 in ["Plots"]:
            for key2 in cfg[key1]:
                title = parse_cfg_plots_title(cfg, key1, key2)
                cfg[key1].rename(key2, title)
                cfg2 = cfg[key1][title]
                for key3 in cfg2:
                    # force keywords to lower case
                    cfg2.rename(key3, key3.lower())
                    cfg3 = cfg2[key3.lower()]
                    cfg3 = parse_cfg_plots_value(key3, cfg3)
                    cfg[key1][title][key3.lower()] = cfg3
        elif key1 in ["Variables"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    if key3 in ["RangeCheck", "DependencyCheck", "DiurnalCheck", "ExcludeDates",
                                "CorrectWindDirection", "ApplyFcStorage",
                                "MergeSeries", "AverageSeries",
                                "LowerCheck", "UpperCheck"]:
                        for key4 in cfg3:
                            # force keywords to lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
                    if key3 in ["ExcludeHours"]:
                        for key4 in cfg3:
                            # force keywords to lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
        else:
            del cfg[key1]
    return cfg

def l2_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: May 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over variables in the control file
    labels_cfg = list(cfg["Variables"].keys())
    for label_cfg in labels_cfg:
        if "DependencyCheck" in cfg["Variables"][label_cfg]:
            source = cfg["Variables"][label_cfg]["DependencyCheck"]["source"]
            vs = pfp_cfg.cfg_string_to_list(source)
            vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
            for rp in renames_pattern:
                llen = len(rp)
                srp = std["rename_pattern"][rp]
                klen = len(srp)
                vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                   (v[:klen] != srp)) else v for v in vs]
            cfg["Variables"][label_cfg]["DependencyCheck"]["source"] = ",".join(vs)
    return cfg

def l2_update_cfg_variable_deprecate(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L2 control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["variables"])
    for label in deprecated:
        if label in cfg["Variables"]:
            cfg["Variables"].pop(label)
    return cfg

def l2_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # loop over the variables in the Variables section of the control file
    for label_cfg in list(cfg["Variables"].keys()):
        if label_cfg in renames_exact:
            new_name = std["rename_exact"][label_cfg]
            cfg["Variables"].rename(label_cfg, new_name)
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg["Variables"].keys()):
            if ((label_cfg[:llen] == rp) and
                (label_cfg[:klen] != srp)):
                new_name = label_cfg.replace(label_cfg[:llen], srp)
                cfg["Variables"].rename(label_cfg, new_name)
    # loop over the variables in the Plots section of the control file
    if "Plots" in list(cfg.keys()):
        plots = list(cfg["Plots"])
        for plot in plots:
            if "variables" in cfg["Plots"][plot]:
                cp = cfg["Plots"][plot]["variables"]
                vs = pfp_cfg.cfg_string_to_list(cp)
                vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                       (v[:klen] != srp)) else v for v in vs]
                cfg["Plots"][plot]["variables"] = ",".join(vs)
            elif "type" in cfg["Plots"][plot]:
                if cfg["Plots"][plot]["type"] == "xy":
                    for axis in ["xseries", "yseries"]:
                        cp = cfg["Plots"][plot][axis]
                        vs = pfp_cfg.cfg_string_to_list(cp)
                        vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                        for rp in renames_pattern:
                            llen = len(rp)
                            srp = std["rename_pattern"][rp]
                            klen = len(srp)
                            vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                               (v[:klen] != srp)) else v for v in vs]
                        cfg["Plots"][plot][axis] = ",".join(vs)
    return cfg

def l3_update_controlfile(cfg):
    """
    Purpose:
     Parse the L3 control file to update the syntax from earlier OFQC/PFP versions
     to the syntax used by this version.
    Usage:
     result = pfp_compliance.update_l3_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L3 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: February 2020
    """
    # get a copy of the original control file object
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        cfg = l3_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L3 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        cfg = l3_update_cfg_options(cfg, std)
        cfg = l3_update_cfg_variable_deprecate(cfg, std)
        cfg = l3_update_cfg_variable_names(cfg, std)
        cfg = l3_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L3 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l3_update_cfg_options(cfg, std):
    """
    Purpose:
     Update an L3 control file Options section.
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: May 2020
    """
    # move items from the General section (if it exists) to the Options section.
    if "General" in cfg:
        if "Options" not in cfg:
            cfg["Options"] = {}
        for item in cfg["General"]:
            cfg["Options"][item] = cfg["General"][item]
        del cfg["General"]
    # update the units in the Options section
    if "Options" in cfg:
        for item in cfg["Options"]:
            if item in ["ApplyFcStorage", "CcUnits", "FcUnits", "ReplaceFcStorage", "DisableFcWPL"]:
                cfg["Options"].rename(item, item.replace("Fc", "Fco2"))
        old_units = list(std["units_map"].keys())
        for item in list(cfg["Options"].keys()):
            if cfg["Options"][item] in old_units:
                cfg["Options"][item] = std["units_map"][cfg["Options"][item]]
    return cfg

def l3_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L3 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: May 2020
    """
    strip_list = ['"', "'", "[", "]"]
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files", "Global", "Soil", "Massman"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            for key2 in cfg[key1]:
                if key2 in ["MassmanCorrection", "UseL2Fluxes",
                            "CO2Units", "Fco2Units", "ApplyFco2Storage",
                            "ReplaceFco2Storage", "DisableFco2WPL",
                            "CcUnits", "FcUnits", "ApplyFcStorage",
                            "ReplaceFcStorage", "DisableFcWPL",
                            "CorrectFgForStorage", "CorrectIndividualFg",
                            "2DCoordRotation", "zms"]:
                    cfg2 = cfg[key1][key2]
                    cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                    cfg[key1][key2] = cfg2
                else:
                    del cfg[key1][key2]
        elif key1 in ["Imports"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    cfg3 = parse_cfg_values(key3, cfg3, strip_list)
                    cfg[key1][key2][key3] = cfg3
        elif key1 in ["Plots"]:
            for key2 in cfg[key1]:
                title = parse_cfg_plots_title(cfg, key1, key2)
                cfg[key1].rename(key2, title)
                cfg2 = cfg[key1][title]
                for key3 in cfg2:
                    # force keywords to lower case
                    cfg2.rename(key3, key3.lower())
                    cfg3 = cfg2[key3.lower()]
                    cfg3 = parse_cfg_plots_value(key3, cfg3)
                    cfg[key1][title][key3.lower()] = cfg3
        elif key1 in ["Variables"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    if key3 in ["RangeCheck", "DependencyCheck", "DiurnalCheck", "ExcludeDates",
                                "ApplyFco2Storage", "MergeSeries", "AverageSeries"]:
                        cfg3 = cfg[key1][key2][key3]
                        for key4 in cfg3:
                            # force keywords to lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
        else:
            del cfg[key1]
    return cfg

def l3_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: May 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over variables in the control file
    labels_cfg = list(cfg["Variables"].keys())
    for label_cfg in labels_cfg:
        for qc in ["DependencyCheck", "MergeSeries", "AverageSeries"]:
            if qc in cfg["Variables"][label_cfg]:
                source = cfg["Variables"][label_cfg][qc]["source"]
                vs = pfp_cfg.cfg_string_to_list(source)
                vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                       (v[:klen] != srp)) else v for v in vs]
                cfg["Variables"][label_cfg][qc]["source"] = ",".join(vs)
            else:
                continue
    return cfg

def l3_update_cfg_variable_deprecate(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L3 control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["variables"])
    for label in deprecated:
        if label in cfg["Variables"]:
            cfg["Variables"].pop(label)
    return cfg

def l3_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # loop over the variables in the Variables section of the control file
    for label_cfg in list(cfg["Variables"].keys()):
        if label_cfg in renames_exact:
            new_name = std["rename_exact"][label_cfg]
            cfg["Variables"].rename(label_cfg, new_name)
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg["Variables"].keys()):
            if ((label_cfg[:llen] == rp) and
                (label_cfg[:klen] != srp)):
                new_name = label_cfg.replace(label_cfg[:llen], srp)
                cfg["Variables"].rename(label_cfg, new_name)
    # loop over the variables in the Plots section of the control file
    if "Plots" in list(cfg.keys()):
        plots = list(cfg["Plots"])
        for plot in plots:
            cp = cfg["Plots"][plot]["variables"]
            vs = pfp_cfg.cfg_string_to_list(cp)
            vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
            for rp in renames_pattern:
                llen = len(rp)
                srp = std["rename_pattern"][rp]
                klen = len(srp)
                vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                   (v[:klen] != srp)) else v for v in vs]
            cfg["Plots"][plot]["variables"] = ",".join(vs)
    return cfg

def l4_update_controlfile(cfg):
    """
    Purpose:
     Parse the L4 control file to make sure the syntax is correct and that the
     control file contains all of the information needed.
    Usage:
     result = pfp_compliance.l4_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L4 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: June 2020
    """
    # get a copy of the original control file object
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        cfg = l4_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L4 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        #cfg = l3_update_cfg_options(cfg, std)
        cfg = l4_update_cfg_variable_deprecate(cfg, std)
        cfg = l4_update_cfg_variable_names(cfg, std)
        cfg = l4_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L4 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l4_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L4 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: June 2020
    """
    strip_list = ['"', "'", "[", "]"]
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files", "Global"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                if key2 in ["MaxGapInterpolate"]:
                    cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                    cfg[key1][key2] = cfg2
                else:
                    del cfg[key1][key2]
        elif key1 in ["Imports"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    cfg3 = parse_cfg_values(key3, cfg3, strip_list)
                    cfg[key1][key2][key3] = cfg3
        elif key1 in ["Drivers"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    if key3 in ["GapFillFromAlternate", "GapFillFromClimatology",
                                "GapFillUsingMDS"]:
                        for key4 in cfg3:
                            cfg4 = cfg[key1][key2][key3][key4]
                            for key5 in cfg4:
                                cfg4.rename(key5, key5.lower())
                                cfg5 = cfg4[key5.lower()]
                                cfg5 = parse_cfg_values(key5, cfg5, strip_list)
                                cfg[key1][key2][key3][key4][key5.lower()] = cfg5
                    elif key3 in ["RangeCheck", "DependencyCheck", "DiurnalCheck",
                                  "ExcludeDates", "MergeSeries"]:
                        # strip out unwanted characters
                        for key4 in cfg3:
                            # force lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg[key1][key2][key3][key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
        elif key1 in ["GUI"]:
            continue
        else:
            del cfg[key1]
    return cfg

def l4_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: June 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over variables in the control file
    labels_cfg = list(cfg["Drivers"].keys())
    for label_cfg in labels_cfg:
        for qc in ["DependencyCheck", "MergeSeries", "AverageSeries"]:
            if qc in cfg["Drivers"][label_cfg]:
                source = cfg["Drivers"][label_cfg][qc]["source"]
                vs = pfp_cfg.cfg_string_to_list(source)
                vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                       (v[:klen] != srp)) else v for v in vs]
                cfg["Drivers"][label_cfg][qc]["source"] = ",".join(vs)
            else:
                continue
        for gfm in ["GapFillFromAlternate", "GapFillFromClimatology", "GapFillUsingMDS"]:
            if gfm in cfg["Drivers"][label_cfg]:
                gfvs = list(cfg["Drivers"][label_cfg][gfm].keys())
                for gfv in gfvs:
                    if gfv in renames_exact:
                        new_name = std["rename_exact"][gfv]
                        cfg["Drivers"][label_cfg][gfm].rename(gfv, new_name)
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    # loop over the variables in the control file
                    for gfv in gfvs:
                        if ((gfv[:llen] == rp) and
                            (gfv[:klen] != srp)):
                            new_name = gfv.replace(gfv[:llen], srp)
                            cfg["Drivers"][label_cfg][gfm].rename(gfv, new_name)
    return cfg

def l4_update_cfg_variable_deprecate(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L4 control file.
    Usage:
    Author: PRI
    Date: June 2020
    """
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["variables"])
    for label in deprecated:
        if label in cfg["Drivers"]:
            cfg["Drivers"].pop(label)
    return cfg

def l4_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: June 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # loop over the variables in the Variables section of the control file
    for label_cfg in list(cfg["Drivers"].keys()):
        if label_cfg in renames_exact:
            new_name = std["rename_exact"][label_cfg]
            cfg["Drivers"].rename(label_cfg, new_name)
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg["Drivers"].keys()):
            if ((label_cfg[:llen] == rp) and
                (label_cfg[:klen] != srp)):
                new_name = label_cfg.replace(label_cfg[:llen], srp)
                cfg["Drivers"].rename(label_cfg, new_name)
    return cfg

def l5_update_controlfile(cfg):
    """
    Purpose:
     Parse the L5 control file to make sure the syntax is correct and that the
     control file contains all of the information needed.
    Usage:
     result = pfp_compliance.l5_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L5 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: June 2020
    """
    # get a copy of the original control file object
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    try:
        cfg = l5_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L5 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        cfg = l5_update_cfg_options(cfg, std)
        cfg = l5_update_cfg_variable_deprecate(cfg, std)
        cfg = l5_update_cfg_variable_names(cfg, std)
        cfg = l5_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L5 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l5_update_cfg_options(cfg, std):
    """
    Purpose:
     Update Options section in L5 control file.
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: June 2020
    """
    # check to see if control file has an Options section
    if "Options" not in cfg:
        return
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over entries in the Options section, exact match
    options = cfg["Options"]
    for option in options:
        opt = cfg["Options"][option]
        if opt in renames_exact:
            cfg["Options"][option] = std["rename_exact"][opt]
    # loop over entries in the Options section, pattern match
    options = cfg["Options"]
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        for option in options:
            opt = str(cfg["Options"][option])
            if ((opt[:llen] == rp) and
                (opt[:klen] != srp)):
                new_opt = opt.replace(opt[:llen], srp)
                cfg["Options"][option] = new_opt
    return cfg

def l5_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L5 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: June 2020
    """
    strip_list = ['"', "'", "[", "]"]
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files", "Global", "ustar_threshold"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            for key2 in cfg[key1]:
                if key2 in ["MaxGapInterpolate", "MaxShortGapLength", "FilterList",
                            "TurbulenceFilter", "DayNightFilter", "AcceptDayTimes",
                            "TruncateToImports"]:
                    cfg2 = cfg[key1][key2]
                    cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                    cfg[key1][key2] = cfg2
                else:
                    del cfg[key1][key2]
        elif key1 in ["Imports"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    cfg3 = parse_cfg_values(key3, cfg3, strip_list)
                    cfg[key1][key2][key3] = cfg3
        elif key1 in ["Fluxes"]:
            # key2 is the variable name
            for key2 in cfg[key1]:
                # key3 is the gap filling method
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    if key3 in ["GapFillUsingSOLO", "GapFillLongSOLO", "GapFillUsingMDS",
                                "GapFillFromClimatology"]:
                        # key4 is the gap fill variable
                        for key4 in cfg3:
                            cfg4 = cfg[key1][key2][key3][key4]
                            for key5 in cfg4:
                                cfg4.rename(key5, key5.lower())
                                cfg5 = cfg4[key5.lower()]
                                cfg5 = parse_cfg_values(key5, cfg5, strip_list)
                                cfg[key1][key2][key3][key4][key5.lower()] = cfg5
                    elif key3 in ["MergeSeries", "RangeCheck", "DependencyCheck",
                                  "DiurnalCheck", "ExcludeDates"]:
                        # strip out unwanted characters
                        for key4 in cfg3:
                            # force lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
        elif key1 in ["GUI"]:
            continue
        else:
            del cfg[key1]
    return cfg

def l5_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: June 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over variables in the control file
    labels_cfg = list(cfg["Fluxes"].keys())
    for label_cfg in labels_cfg:
        for qc in ["DependencyCheck", "MergeSeries", "AverageSeries"]:
            if qc in cfg["Fluxes"][label_cfg]:
                source = cfg["Fluxes"][label_cfg][qc]["source"]
                vs = pfp_cfg.cfg_string_to_list(source)
                vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                       (v[:klen] != srp)) else v for v in vs]
                cfg["Fluxes"][label_cfg][qc]["source"] = ",".join(vs)
            else:
                continue
        for gfm in ["GapFillUsingSOLO", "GapFillLongSOLO", "GapFillUsingMDS"]:
            if gfm in cfg["Fluxes"][label_cfg]:
                # loop over the variables for this gap filling method and rename if necessary
                gfvs = list(cfg["Fluxes"][label_cfg][gfm].keys())
                for gfv in gfvs:
                    if gfv in renames_exact:
                        new_name = std["rename_exact"][gfv]
                        cfg["Fluxes"][label_cfg][gfm].rename(gfv, new_name)

                # now loop over the drivers and rename if necessary
                gfvs = list(cfg["Fluxes"][label_cfg][gfm].keys())
                for gfv in gfvs:
                    drivers = cfg["Fluxes"][label_cfg][gfm][gfv]["drivers"]
                    for re in renames_exact:
                        if re in drivers:
                            drivers = drivers.replace(re, std["rename_exact"][re])
                    cfg["Fluxes"][label_cfg][gfm][gfv]["drivers"] = drivers

                # rename if first characters of variable name match pattern
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    # loop over the variables in the control file
                    gfvs = list(cfg["Fluxes"][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        if ((gfv[:llen] == rp) and
                            (gfv[:klen] != srp)):
                            new_name = gfv.replace(gfv[:llen], srp)
                            cfg["Fluxes"][label_cfg][gfm].rename(gfv, new_name)
                # loop over the drivers
                gfvs = list(cfg["Fluxes"][label_cfg][gfm].keys())
                for gfv in gfvs:
                    drivers = cfg["Fluxes"][label_cfg][gfm][gfv]["drivers"]
                    for rp in renames_pattern:
                        if rp in drivers:
                            drivers = drivers.replace(rp, std["rename_pattern"][rp])
    return cfg

def l5_update_cfg_variable_deprecate(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L5 control file.
    Usage:
    Author: PRI
    Date: June 2020
    """
    deprecated = pfp_cfg.cfg_string_to_list(std["deprecated"]["variables"])
    for label in deprecated:
        if label in cfg["Fluxes"]:
            cfg["Fluxes"].pop(label)
    return cfg

def l5_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: June 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["rename_exact"].keys())
    # loop over the variables in the Fluxes section of the control file
    for label_cfg in list(cfg["Fluxes"].keys()):
        if label_cfg in renames_exact:
            new_name = std["rename_exact"][label_cfg]
            cfg["Fluxes"].rename(label_cfg, new_name)
    # rename pattern matches
    renames_pattern = list(std["rename_pattern"].keys())
    for rp in renames_pattern:
        llen = len(rp)
        srp = std["rename_pattern"][rp]
        klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg["Fluxes"].keys()):
            if ((label_cfg[:llen] == rp) and
                (label_cfg[:klen] != srp)):
                new_name = label_cfg.replace(label_cfg[:llen], srp)
                cfg["Fluxes"].rename(label_cfg, new_name)
    return cfg

def l6_update_controlfile(cfg):
    """
    Purpose:
     Parse the L6 control file to make sure the syntax is correct and that the
     control file contains all of the information needed.
    Usage:
     result = pfp_compliance.l6_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the L6 control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: June 2020
    """
    # get a copy of the original control file object
    cfg_original = copy.deepcopy(cfg)
    # initialise the return logical
    ok = True
    # check to see if we have an old style L6 control file
    if "ER" in list(cfg.keys()) or "NEE" in list(cfg.keys()) or "GPP" in list(cfg.keys()):
        ok = False
        msg = "This is an old version of the L6 control file.\n"
        msg = msg + "Close the L6 control file and create a new one from\n"
        msg = msg + "the template in PyFluxPro/controlfiles/template/L6."
        result = pfp_gui.MsgBox_Quit(msg, title="Critical")
        return ok
    try:
        cfg = l6_update_cfg_syntax(cfg)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L6 control file syntax"
    # check to see if we can load the nc_cleanup.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "cfg_update.txt")
        std = pfp_io.get_controlfilecontents(stdname)
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # clean up the variable names
    try:
        #cfg = l5_update_cfg_options(cfg, std)
        #cfg = l5_update_cfg_variable_deprecate(cfg, std)
        cfg = l6_update_cfg_variable_names(cfg, std)
        cfg = l6_update_cfg_variable_attributes(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L6 control file contents"
    if ok:
        # check to see if the control file object has been changed
        if cfg != cfg_original:
            # and save it if it has changed
            file_name = os.path.basename(cfg.filename)
            logger.info(" Updated and saved control file " + file_name)
            cfg.write()
    else:
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return ok

def l6_update_cfg_syntax(cfg):
    """
    Purpose:
     Update an L6 control file from the old-style (pre PFP V1.0) syntax to the
     new-style syntax (post PFP V1.0).
    Usage:
    Side effects:
     Returns a modified control file object.
    Author: PRI
    Date: June 2020
    """
    strip_list = ['"', "'", "[", "]"]
    for key1 in cfg:
        if key1 in ["level"]:
            cfg[key1] = parse_cfg_values(key1, cfg[key1], strip_list)
        elif key1 in ["Files", "Global"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            # no options section in L6 control file yet
            pass
        elif key1 in ["EcosystemRespiration"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                for key3 in cfg2:
                    cfg3 = cfg[key1][key2][key3]
                    if key3 in ["ERUsingSOLO", "ERUsingLloydTaylor", "ERUsingLasslop"]:
                        for key4 in cfg3:
                            cfg4 = cfg[key1][key2][key3][key4]
                            for key5 in cfg4:
                                cfg4.rename(key5, key5.lower())
                                cfg5 = cfg4[key5.lower()]
                                cfg5 = parse_cfg_values(key5, cfg5, strip_list)
                                cfg[key1][key2][key3][key4][key5.lower()] = cfg5
                    elif key3 in ["MergeSeries"]:
                        # strip out unwanted characters
                        for key4 in cfg[key1][key2][key3]:
                            # force lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
        elif key1 in ["NetEcosystemExchange", "GrossPrimaryProductivity"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    cfg3 = parse_cfg_values(key3, cfg3, strip_list)
                    cfg[key1][key2][key3] = cfg3
        elif key1 in ["GUI"]:
            continue
        else:
            del cfg[key1]
    return cfg

def l6_update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: June 2020
    """
    # exact and pattern rename lists
    renames_exact = list(std["rename_exact"].keys())
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over the variables in the EcosystemRespiration section of the control file
    for method in ["EcosystemRespiration"]:
        labels_cfg = list(cfg[method].keys())
        for label_cfg in labels_cfg:
            for qc in ["MergeSeries", "AverageSeries"]:
                if qc in cfg[method][label_cfg]:
                    source = cfg[method][label_cfg][qc]["source"]
                    vs = pfp_cfg.cfg_string_to_list(source)
                    vs = [std["rename_exact"][v] if v in renames_exact else v for v in vs]
                    for rp in renames_pattern:
                        llen = len(rp)
                        srp = std["rename_pattern"][rp]
                        klen = len(srp)
                        vs = [v.replace(v[:llen], srp) if ((v[:llen] == rp) and
                                                           (v[:klen] != srp)) else v for v in vs]
                    cfg[method][label_cfg][qc]["source"] = ",".join(vs)
                else:
                    continue
            for gfm in ["ERUsingSOLO", "ERUsingLloydTaylor", "ERUsingLasslop"]:
                if gfm in cfg[method][label_cfg]:
                    # loop over the variables for this gap filling method and rename if necessary
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        if gfv in renames_exact:
                            new_name = std["rename_exact"][gfv]
                            cfg[method][label_cfg][gfm].rename(gfv, new_name)

                    # now loop over the drivers and rename if necessary
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        item = cfg[method][label_cfg][gfm][gfv]["drivers"]
                        drivers = pfp_cfg.cfg_string_to_list(item)
                        for driver in list(drivers):
                            for re in renames_exact:
                                if driver == re:
                                    new_name = std["rename_exact"][re]
                                    drivers[drivers.index(driver)] = new_name
                        cfg[method][label_cfg][gfm][gfv]["drivers"] = ",".join(drivers)

                    # rename if first characters of variable name match pattern
                    for rp in renames_pattern:
                        llen = len(rp)
                        srp = std["rename_pattern"][rp]
                        klen = len(srp)
                        # loop over the variables in the control file
                        gfvs = list(cfg[method][label_cfg][gfm].keys())
                        for gfv in gfvs:
                            if ((gfv[:llen] == rp) and
                                (gfv[:klen] != srp)):
                                new_name = gfv.replace(gfv[:llen], srp)
                                cfg[method][label_cfg][gfm].rename(gfv, new_name)

                    # loop over the drivers
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        item = cfg[method][label_cfg][gfm][gfv]["drivers"]
                        drivers = pfp_cfg.cfg_string_to_list(item)
                        for driver in list(drivers):
                            for rp in renames_pattern:
                                llen = len(rp)
                                srp = std["rename_pattern"][rp]
                                klen = len(srp)
                                if ((driver[:llen] == rp) and
                                    (driver[:klen] != srp)):
                                    new_name = opt.replace(driver[:llen], srp)
                                    drivers[drivers.index(driver)] = new_name
                        cfg[method][label_cfg][gfm][gfv]["drivers"] = ",".join(drivers)


    # loop over the variables in the NEE and GPP sections
    for method in ["NetEcosystemExchange", "GrossPrimaryProductivity"]:
        labels_cfg = list(cfg[method].keys())
        for label_cfg in labels_cfg:
            # first loop is over the keys
            variables = list(cfg[method][label_cfg].keys())
            for var in variables:
                opt = cfg[method][label_cfg][var]
                # exact matches
                if var in renames_exact:
                    new_name = std["rename_exact"][var]
                    cfg[method][label_cfg].rename(var, new_name)
                # pattern matches
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    if ((var[:llen] == rp) and
                        (var[:klen] != srp)):
                        new_name = var.replace(var[:llen], srp)
                        cfg[method][label_cfg].rename(var, new_name)
            # second loop is over the values in case any of the keys were renamed
            variables = list(cfg[method][label_cfg].keys())
            for var in variables:
                opt = cfg[method][label_cfg][var]
                # exact matches
                if opt in renames_exact:
                    new_name = std["rename_exact"][opt]
                    cfg[method][label_cfg][var] = new_name
                # pattern matches
                for rp in renames_pattern:
                    llen = len(rp)
                    srp = std["rename_pattern"][rp]
                    klen = len(srp)
                    if ((opt[:llen] == rp) and
                        (opt[:klen] != srp)):
                        new_name = opt.replace(opt[:llen], srp)
                        cfg[method][label_cfg][var] = new_name
    return cfg

def l6_update_cfg_variable_names(cfg, std):
    """
    Purpose:
     Update the variable names according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: June 2020
    """
    # exact and pattern rename lists
    renames_exact = list(std["rename_exact"].keys())
    renames_pattern = list(std["rename_pattern"].keys())
    # loop over the variables in the EcosystemRespiration, NetEcosystemExchange and
    # GrossPrimaryProductivity section of the control file
    for method in ["EcosystemRespiration", "NetEcosystemExchange", "GrossPrimaryProductivity"]:
        # rename exact matches
        for label_cfg in list(cfg[method].keys()):
            if label_cfg in renames_exact:
                new_name = std["rename_exact"][label_cfg]
                cfg[method].rename(label_cfg, new_name)
        # rename pattern matches
        for rp in renames_pattern:
            llen = len(rp)
            srp = std["rename_pattern"][rp]
            klen = len(srp)
            # loop over the variables in the control file
            for label_cfg in list(cfg[method].keys()):
                if ((label_cfg[:llen] == rp) and
                    (label_cfg[:klen] != srp)):
                    new_name = label_cfg.replace(label_cfg[:llen], srp)
                    cfg[method].rename(label_cfg, new_name)
    return cfg

def parse_cfg_plots_title(cfg, key1, key2):
    """ Parse the [Plots] section for a title."""
    title = key2
    for item in list(cfg[key1][key2].keys()):
        if item.lower() == "title":
            title = cfg[key1][key2][item]
            del cfg[key1][key2][item]
            break
    strip_list = ['"', "'"]
    for c in strip_list:
        if c in title:
            title = title.replace(c, "")
    return title

def parse_cfg_plots_value(k, v):
    """ Parse the [Plots] section keys to remove unnecessary characters."""
    if k.lower() in ["variables", "type", "xseries", "yseries"]:
        if ("[" in v) and ("]" in v):
            v = v.replace("[", "").replace("]", "")
    strip_list = [" ", '"', "'"]
    for c in strip_list:
        if c in v:
            v = v.replace(c, "")
    return v

def parse_cfg_values(k, v, strip_list):
    """ Parse key values to remove unnecessary characters."""
    # strip unwanted characters
    for c in strip_list:
        if c in v:
            v = v.replace(c, "")
    if k in ["file_path", "plot_path"] and "browse" not in v:
        if os.path.join(str(v), "") != v:
            v = os.path.join(str(v), "")
    return v

def parse_cfg_variables_excludehours(k, v):
    """ Parse value from ExcludeHours to remove unwanted characters"""
    v = v.replace("[", "").replace("]", "")
    if k in ["ExcludeHours"]:
        strip_list = ["'", '"']
    # strip unwanted characters
    for c in strip_list:
        if c in v:
            v = v.replace(c, "")
    return v

def parse_cfg_variables_value(k, v):
    """ Parse value from control file to remove unnecessary characters."""
    try:
        # check to see if it is a number
        r = float(v)
    except ValueError as e:
        if ("[" in v) and ("]" in v) and ("*" in v):
            # old style of [value]*12
            v = v[v.index("[")+1:v.index("]")]
        elif ("[" in v) and ("]" in v) and ("*" not in v):
            # old style of [1,2,3,4,5,6,7,8,9,10,11,12]
            v = v.replace("[", "").replace("]", "")
    # remove white space and quotes
    if k in ["ExcludeDates", "CorrectWindDirection", "LowerCheck", "UpperCheck"]:
        # don't remove white space between date and time
        strip_list = ['"', "'"]
    elif k in ["Attr"]:
        # don't remove white space from variable attributes
        strip_list = ['"', "'"]
    else:
        strip_list = [" ", '"', "'"]
    # strip unwanted characters
    for c in strip_list:
        if c in v:
            v = v.replace(c, "")
    return v

def nc_update(cfg):
    """
    Purpose:
     Update an OFQC-style netCDF by changing variable names and attributes.
    Usage:
    Author: PRI
    Date: October 2018
    """
    # get the input file path
    nc_file_path = pfp_io.get_infilenamefromcf(cfg)
    msg = " Converting file " + os.path.split(nc_file_path)[1]
    logger.info(msg)
    # read the input file
    ds1 = pfp_io.nc_read_series(nc_file_path)
    if ds1.returncodes["value"] != 0: return
    # update the variable names
    change_variable_names(cfg, ds1)
    # make sure there are Ws and Wd series
    copy_ws_wd(ds1)
    # make sure we have all the variables we want ...
    ds2 = include_variables(cfg, ds1)
    # ... but not the ones we don't
    exclude_variables(cfg, ds2)
    # update the variable units
    change_variable_units(cfg, ds2)
    # make sure there are Ws and Wd series
    copy_ws_wd(ds2)
    ## make sure we have all the variables we want ...
    #ds2 = include_variables(cfg, ds1)
    ## ... but not the ones we don't
    #exclude_variables(cfg, ds2)
    # update the global attributes
    change_global_attributes(cfg, ds2)
    # update the variable attributes
    change_variable_attributes(cfg, ds2)
    # Fc single point storage
    consistent_Fco2_storage(ds2, os.path.split(nc_file_path)[1])
    # rename the original file to prevent it being overwritten
    t = time.localtime()
    rundatetime = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5]).strftime("%Y%m%d%H%M")
    new_ext = "_" + rundatetime + ".nc"
    # add the current local datetime the base file name
    new_file_path = nc_file_path.replace(".nc", new_ext)
    msg = " Renaming original file to " + os.path.split(new_file_path)[1]
    logger.info(msg)
    # ... and rename the base file to preserve it
    os.rename(nc_file_path, new_file_path)
    # write the updated file
    nc_file = pfp_io.nc_open_write(nc_file_path)
    if nc_file is None: return
    pfp_io.nc_write_series(nc_file, ds2)

    return 0

