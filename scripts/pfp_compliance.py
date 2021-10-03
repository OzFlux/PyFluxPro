# standard modules
import copy
import inspect
import logging
import os
import platform
import subprocess
import traceback
# 3rd party modules
# PFP modules
from scripts import pfp_func_units
from scripts import pfp_func_stats
from scripts import pfp_gui
from scripts import pfp_io
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def CheckCFCompliance(nc_file_uri):
    """
    Purpose:
     Run the CF Conventions checker on a netCDF file.
    Side effects:
     Creates 3 text files in the same folder as the netCDF file
     and with the following suffixes:
      _cfchecker.all - contains all messages from cfchecker
      _cfchecker.errors - contains all ERROR messages from cfchecker
      _cfchecker.warnings - contains all WARNING messages from cfchecker
    Author: PRI
    Date: July 2021
    """
    msg = " Checking CF compliance for " + os.path.basename(nc_file_uri)
    logger.info(msg)
    try:
        # cfchecker log file for all messages
        cf_file_uri = nc_file_uri.replace(".nc", "_cfchecker.all")
        # cfchecker log file for error messages
        cf_error_uri = cf_file_uri.replace("_cfchecker.all", "_cfchecker.errors")
        # cfchecker log file for warning messages
        cf_warning_uri = cf_file_uri.replace("_cfchecker.all", "_cfchecker.warnings")
        # with a little work, cfchecker could be used directly in PFP
        cmd = ["cfchecks", "-v 1.8", nc_file_uri]
        with open(cf_file_uri, "w") as f:
            # run cfchecks as a subprocess at this stage
            subprocess.run(cmd, stdout=f)
        # parse the output of cfchecker and write warnings and errors to separate files
        # error file
        error_file = open(cf_error_uri, "w")
        n_errors = 0
        # warning file
        warning_file = open(cf_warning_uri, "w")
        n_warnings = 0
        # read the cfchecker log file
        with open(cf_file_uri) as f:
            for line in f:
                if "CHECKING NetCDF FILE" in line:
                    # line containing the file name
                    error_file.write(line)
                    warning_file.write(line)
                    continue
                if "Checking variable:" in line:
                    # line containing the variable name
                    parts = line.split(" ")
                    variable_name = parts[2].rstrip()
                    continue
                if (("ERROR:" in line) and ('variable_name' in locals())):
                    # line containing error messages
                    line_out = " " + variable_name + " " + line
                    error_file.write(line_out)
                    n_errors += 1
                    continue
                if (("WARN:" in line) and ('variable_name' in locals())):
                    # line containing warning messages
                    line_out = " " + variable_name + " " + line
                    warning_file.write(line_out)
                    n_warnings += 1
                    continue
        # close the error and warning files
        error_file.close()
        warning_file.close()
        # tell the user what we found
        if n_errors > 0:
            # number of errors found by cfchecker
            msg = " CF checks returned " + str(n_errors) + " errors"
            logger.error(msg)
            msg = "  Check " + cf_error_uri + " for details"
            logger.error(msg)
        if n_warnings > 0:
            # number of warnings found by cfchecker
            msg = " CF checks returned " + str(n_warnings) + " warnings"
            logger.warning(msg)
            msg = "  Check " + cf_warning_uri + " for details"
            logger.warning(msg)
        # info message to user about errors and warnings
        msg = " CF checks returned " + str(n_errors) + " errors and "
        msg += str(n_warnings) + " warnings"
        logger.info(msg)
    except Exception:
        msg = "Error during CF compliance check of " + os.path.split(nc_file_uri)[1]
        logger.error(msg)
        error_message = traceback.format_exc()
        logger.error(error_message)
    return

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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the climatology control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_rename(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the climatology control file contents"
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    # force the level to "cpd_barr", early versions used "cpd1"
    cfg["level"] = "cpd_barr"
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the CPD (Barr) control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_rename(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the CPD (Barr) control file contents"
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    # force the level to "cpd_mchugh", early versions used "cpd2"
    cfg["level"] = "cpd_mchugh"
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the CPD (McHugh) control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_rename(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the CPD (McHugh) control file contents"
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

def cpd_mcnew_update_controlfile(cfg):
    """
    Purpose:
     Parse the CPD (McNew) control file to update the syntax from earlier OFQC/PFP
     versions to the syntax used by this version.
    Usage:
     result = pfp_compliance.cpd_mcnew_update_controlfile(cfg)
     where cfg is a ConfigObj object
           result is True if the concatenate control file was updated successfully
                     False if it couldn't be updated
    Side effects:
    Author: PRI
    Date: May 2021
    """
    # copy the control file
    cfg_original = copy.deepcopy(cfg)
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    # force the level to "cpd_mcnew", early versions used "cpd3"
    cfg["level"] = "cpd_mcnew"
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the CPD (McNew) control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_rename(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the CPD (McNew) control file contents"
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
    # force the level to "mpt"
    cfg["level"] = "mpt"
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the MPT control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_rename(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the MPT control file contents"
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the concatenation control file syntax"
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
            msg = " File not found (" + os.path.basename(file_name) + ")"
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
    inc["time_coverage_start"] = []
    inc["time_coverage_end"] = []
    inc["chrono_files"] = []
    inc["labels"] = []
    inc["attributes"] = ["height", "instrument", "long_name", "standard_name",
                         "statistic_type", "units", "valid_range"]
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
    l1ire["Files"]["in_headerrow"] = int(cf["Files"]["in_headerrow"])
    l1ire["Files"]["in_firstdatarow"] = int(cf["Files"]["in_firstdatarow"])
    # get the global attributes
    l1ire["Global"] = copy.deepcopy(cf["Global"])
    # get the variables
    l1ire["Variables"] = copy.deepcopy(cf["Variables"])
    # checks of the 'Variables' sections would go here
    for label in list(l1ire["Variables"].keys()):

        # check the 'xl' subsection is present
        ok = True
        if (("xl" not in list(l1ire["Variables"][label].keys())) and
            ("csv" not in list(l1ire["Variables"][label].keys()))):
            # the xl or csv subsections can be missing if Function is present
            if "Function" in list(l1ire["Variables"][label].keys()):
                pass
            else:
                msg = " Skipping " + label + " (subsections 'xl', 'csv' or 'Function' not found)"
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
        elif "csv" in list(l1ire["Variables"][label].keys()):
            for item in ["name"]:
                if item not in list(l1ire["Variables"][label]["csv"].keys()):
                    msg = " Skipping " + label + " ('" + item + "' not found in 'csv' subsection)"
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
        for item in ["long_name"]:
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
        ds.returncodes["value"] = l3_info["status"]["value"]
        ds.returncodes["message"] = l3_info["status"]["message"]
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
    in_file_name = pfp_io.get_infilenamefromcf(cfg)
    if not os.path.isfile(in_file_name):
        msg = " Input file " + os.path.basename(in_file_name) + " not found"
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
    # *******************
    # *** check units ***
    # *******************
    # the lists of units below need to be abstracted from the code and reconciled with
    # the contents of controlfiles/standard/update_control_files.txt
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
             "misc": ["1", "V", "none", "ohm"]}
    ok_units = []
    for key in list(units.keys()):
        ok_units += units[key]
    ok_units = list(set(ok_units))
    for label in list(cfg["Variables"].keys()):
        if "units" not in cfg["Variables"][label]["Attr"]:
            continue
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard",
                               "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred updating the L1 control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_global_attributes(cfg, std)
        cfg = update_cfg_variables_deprecated(cfg, std)
        cfg = update_cfg_variables_rename(cfg, std)
        cfg = l1_update_cfg_variables_attributes(cfg, std)
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

def update_cfg_global_attributes(cfg, std):
    """
    Purpose:
     Update the global attributes according to the rules in the standard control file.
    Usage:
    Author: PRI
    Date: May 2020
          June 2021 rewrite for DSA compliance, phase 1
    """
    stdga = std["Global_attributes"]
    # check for the essential global attributes
    essentials = pfp_utils.string_to_list(stdga["essentials"]["global"])
    for gattr in essentials:
        if gattr not in cfg["Global"]:
            cfg["Global"][gattr] = ""
    # remove deprecated global attributes
    deprecated = pfp_utils.string_to_list(stdga["deprecated"]["global"])
    for gattr in deprecated:
        if gattr in cfg["Global"]:
            cfg["Global"].pop(gattr)
    # rename global attributes
    rename_exact = list(stdga["rename_exact"].keys())
    for old_name in rename_exact:
        if old_name in cfg["Global"]:
            new_name = stdga["rename_exact"][old_name]
            cfg["Global"][new_name] = cfg["Global"].pop(old_name)
    # add or change global attributes as required
    force = list(stdga["force"].keys())
    for item in force:
        cfg["Global"][item] = stdga["force"][item]
    return cfg

def update_cfg_syntax(cfg, std):
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
    for key1 in cfg:
        if key1 in ["level", "GUI"]:
            continue
        elif key1 in ["Files", "Global", "Soil", "Massman", "ustar_threshold"]:
            for key2 in cfg[key1]:
                cfg2 = cfg[key1][key2]
                cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                cfg[key1][key2] = cfg2
        elif key1 in ["Options"]:
            if cfg["level"] in list(std["Options"].keys()):
                options = str(std["Options"][cfg["level"]])
                options = pfp_utils.string_to_list(options.replace("\n", ""))
                for key2 in cfg[key1]:
                    if key2 in options:
                        cfg2 = cfg[key1][key2]
                        cfg2 = parse_cfg_values(key2, cfg2, strip_list)
                        cfg[key1][key2] = cfg2
                    else:
                        del cfg[key1][key2]
            else:
                msg = " No entry for level " + cfg["level"] + " in std control file"
                logger.warning(msg)
        elif key1 in ["Imports"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    cfg3 = parse_cfg_values(key3, cfg3, strip_list)
                    cfg[key1][key2][key3] = cfg3
        elif key1 in ["Variables", "Drivers", "Fluxes"]:
            for key2 in cfg[key1]:
                for key3 in cfg[key1][key2]:
                    cfg3 = cfg[key1][key2][key3]
                    if key3 in ["xl", "csv", "Attr"]:
                        # L1 control file
                        for key4 in cfg3:
                            # for keywords to lower case
                            if key4.lower() != key4:
                                cfg3[key4.lower()] = cfg3.pop(key4)
                            cfg[key1][key2][key3][key4.lower()] = parse_cfg_variables_value(key3, cfg3[key4.lower()])
                    if key3 in ["RangeCheck", "DependencyCheck", "DiurnalCheck", "ExcludeDates",
                                "CorrectWindDirection", "ApplyFco2Storage",
                                "MergeSeries", "AverageSeries",
                                "LowerCheck", "UpperCheck"]:
                        # L2, L3, L4 and L5 control files
                        for key4 in cfg3:
                            # force keywords to lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
                    if key3 in ["ExcludeHours"]:
                        # L2, L3, L4 and L5 control files
                        for key4 in cfg3:
                            # force keywords to lower case
                            cfg3.rename(key4, key4.lower())
                            cfg4 = cfg3[key4.lower()]
                            cfg4 = parse_cfg_variables_value(key3, cfg4)
                            cfg[key1][key2][key3][key4.lower()] = cfg4
                    if key3 in ["GapFillFromAlternate", "GapFillFromClimatology",
                                "GapFillUsingMDS"]:
                        # L4 control file
                        for key4 in cfg3:
                            cfg4 = cfg3[key4]
                            for key5 in cfg4:
                                cfg4.rename(key5, key5.lower())
                                cfg5 = cfg4[key5.lower()]
                                cfg5 = parse_cfg_values(key5, cfg5, strip_list)
                                cfg[key1][key2][key3][key4][key5.lower()] = cfg5
                    if key3 in ["GapFillUsingSOLO", "GapFillLongSOLO", "GapFillUsingMDS",
                                "GapFillFromClimatology"]:
                        # L5 control file
                        for key4 in cfg3:
                            cfg4 = cfg[key1][key2][key3][key4]
                            for key5 in cfg4:
                                cfg4.rename(key5, key5.lower())
                                cfg5 = cfg4[key5.lower()]
                                cfg5 = parse_cfg_values(key5, cfg5, strip_list)
                                cfg[key1][key2][key3][key4][key5.lower()] = cfg5
                    if cfg["level"] in ["climatology", "cpd_barr", "cpd_mchugh"]:
                        # climatology, CPD (Barr), CPD (McHugh) control file
                        cfg[key1][key2][key3] = parse_cfg_values(key3, cfg3, strip_list)
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
        else:
            del cfg[key1]
    return cfg

def l1_update_cfg_variables_attributes(cfg, std):
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
    vattrs_essential = pfp_utils.string_to_list(std["Variables"]["essential"]["attributes"])
    stdvd = std["Variables"]["deprecated"]
    attributes_deprecated = pfp_utils.string_to_list(stdvd["attributes"])
    units_deprecated = pfp_utils.string_to_list(stdvd["units"])
    height_deprecated = pfp_utils.string_to_list(stdvd["height"])
    miscellaneous_deprecated = pfp_utils.string_to_list(stdvd["miscellaneous"])
    standard_name_deprecated = pfp_utils.string_to_list(stdvd["standard_name"])
    # list of standard attribute values
    labels_std = list(std["Variables"]["attributes"].keys())
    labels_cfg = list(cfg["Variables"].keys())
    # add any essential variable attributes that are missing, deprecate those no longer used
    for label in labels_cfg:
        vattrs_cfg = list(cfg["Variables"][label]["Attr"].keys())
        for vattr in vattrs_essential:
            if vattr not in vattrs_cfg:
                cfg["Variables"][label]["Attr"][vattr] = label
        vattrs_cfg = list(cfg["Variables"][label]["Attr"].keys())
        for vattr in attributes_deprecated:
            if vattr in vattrs_cfg:
                del cfg["Variables"][label]["Attr"][vattr]
        vattrs_cfg = list(cfg["Variables"][label]["Attr"].keys())
        for vattr in vattrs_cfg:
            value = cfg["Variables"][label]["Attr"][vattr]
            if vattr == "units":
                if ((value in units_deprecated) or
                    (len(str(value)) == 0)):
                    cfg["Variables"][label]["Attr"].pop(vattr)
            elif vattr == "height":
                if ((value in height_deprecated) or
                    (len(str(value)) == 0)):
                    cfg["Variables"][label]["Attr"].pop(vattr)
            elif vattr == "standard_name":
                if ((value in standard_name_deprecated) or
                    (len(str(value)) == 0)):
                    cfg["Variables"][label]["Attr"].pop(vattr)
            elif len(str(value)) == 0:
                cfg["Variables"][label]["Attr"].pop(vattr)
            else:
                if ((value in miscellaneous_deprecated) or
                    (len(str(value)) == 0)):
                    cfg["Variables"][label]["Attr"].pop(vattr)
    # coerce units into a standard form
    old_units = list(std["Variables"]["units_map"].keys())
    new_units = [std["Variables"]["units_map"][o] for o in old_units]
    ok_units = list(set(old_units + new_units))
    for label_cfg in labels_cfg:
        if "units" not in cfg["Variables"][label_cfg]["Attr"]:
            continue
        cfg_units = cfg["Variables"][label_cfg]["Attr"]["units"]
        if cfg_units not in ok_units:
            msg = " Unrecognised units " + cfg_units + " for variable " + label_cfg
            logger.warning(msg)
            continue
        if cfg_units.lower() == "none":
            cfg["Variables"][label_cfg]["Attr"].pop("units")
        elif len(cfg_units) == 0:
            cfg["Variables"][label_cfg]["Attr"].pop("units")
        elif cfg_units in old_units:
            cfg["Variables"][label_cfg]["Attr"]["units"] = std["Variables"]["units_map"][cfg_units]
        else:
            pass
    # force some variable attributes to particular values
    for label_std in labels_std:
        # length of the label stub in the standard control file
        llen = len(label_std)
        # pointer to attributes in standard control file
        attr_std = std["Variables"]["attributes"][label_std]
        # list of permitted units for variables that match this stub
        units_std = pfp_utils.string_to_list(attr_std["units"])
        # loop over variables in the L1 control file
        labels_cfg = [l for l in list(cfg["Variables"].keys()) if l[:len(label_std)] == label_std]
        for label_cfg in labels_cfg:
            # pointer to attributes in user control file
            attr_cfg = cfg["Variables"][label_cfg]["Attr"]
            # units string given in the L1 control file
            if "units" in attr_cfg:
                units_cfg = attr_cfg["units"]
            else:
                units_cfg = None
            if ((label_cfg[:llen] == label_std) and (units_cfg in units_std) and
                (label_cfg[-3:] != "_Sd") and (label_cfg[-3:] != "_Vr")):
                # the first letters and units match so update the long_name
                attr_cfg["long_name"] = attr_std["long_name"]
                # update the statistic type
                attr_cfg["statistic_type"] = attr_std["statistic_type"]
                if ("standard_name" in attr_std):
                    # list of standard names for variables that match this stub
                    names_std = pfp_utils.string_to_list(attr_std["standard_name"])
                    # update the standard name
                    idx = units_std.index(units_cfg)
                    if idx > len(names_std):
                        idx = 0
                    attr_cfg["standard_name"] = names_std[idx]
            if (label_cfg[:llen] == label_std) and (label_cfg[-3:] == "_Sd"):
                attr_cfg["long_name"] = attr_std["long_name"]
                attr_cfg["statistic_type"] = "standard deviation"
            if (label_cfg[:llen] == label_std) and (label_cfg[-3:] == "_Vr"):
                attr_cfg["long_name"] = attr_std["long_name"]
                attr_cfg["statistic_type"] = "variance"
    ## remove undefined attributes
    #labels_cfg = list(cfg["Variables"].keys())
    #for label_cfg in labels_cfg:
        #for attr in ["units", "standard_name"]:
            #if attr in cfg["Variables"][label_cfg]["Attr"]:
                #if cfg["Variables"][label_cfg]["Attr"][attr] in vattrs_deprecated_values:
                    #cfg["Variables"][label_cfg]["Attr"].pop(attr)
    return cfg

def update_cfg_variables_deprecated(cfg, std):
    """
    Purpose:
     Remove deprecated variables from L1 control file.
    Usage:
    Author: PRI
    Date: May 2020
    """
    stdvd = std["Variables"]["deprecated"]
    deprecated = pfp_utils.string_to_list(stdvd["variables"])
    for item in ["Variables", "Drivers", "Fluxes"]:
        if item in cfg:
            labels = sorted(list(cfg[item].keys()))
            section_name = item
    for label in deprecated:
        if label in labels:
            cfg[section_name].pop(label)
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L2 control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_deprecated(cfg, std)
        cfg = update_cfg_variables_rename(cfg, std)
        cfg = update_cfg_variable_attributes(cfg, std)
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

def update_cfg_variable_attributes(cfg, std):
    """
    Purpose:
     Update the variable attributes according to the rules in the standard control file.
      - rename variables in DependencyCheck
    Usage:
    Author: PRI
    Date: May 2020
    """
    # rename exact variable name matches
    renames_exact = list(std["Variables"]["rename_exact"].keys())
    # rename pattern matches
    renames_pattern = list(std["Variables"]["rename_pattern"].keys())
    # find out the section name
    for item in ["Variables", "Drivers", "Fluxes"]:
        if item in cfg:
            section_name = item
    # loop over variables in the control file
    labels_cfg = list(cfg[section_name].keys())
    for label_cfg in labels_cfg:
        for qc in ["DependencyCheck", "MergeSeries", "AverageSeries"]:
            if qc in cfg[section_name][label_cfg]:
                source = cfg[section_name][label_cfg][qc]["source"]
                vs = pfp_utils.string_to_list(source)
                vs = [std["Variables"]["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    vs = [v.replace(v.split("_")[0], srp)
                          if (v.split("_")[0] == rp)
                          else v for v in vs]
                cfg[section_name][label_cfg][qc]["source"] = ",".join(vs)
            else:
                continue

        for gfm in ["GapFillFromAlternate", "GapFillFromClimatology"]:
            if gfm in cfg[section_name][label_cfg]:
                gfvs = list(cfg[section_name][label_cfg][gfm].keys())
                for gfv in gfvs:
                    if gfv in renames_exact:
                        new_name = std["Variables"]["rename_exact"][gfv]
                        cfg[section_name][label_cfg][gfm].rename(gfv, new_name)
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    # loop over the variables in the control file
                    for gfv in gfvs:
                        if (gfv.split("_")[0] == rp):
                            new_name = gfv.replace(gfv.split("_")[0], srp)
                            cfg[section_name][label_cfg][gfm].rename(gfv, new_name)
            else:
                continue

        for gfm in ["GapFillUsingSOLO", "GapFillLongSOLO", "GapFillUsingMDS"]:
            if gfm in cfg[section_name][label_cfg]:
                # loop over the variables for this gap filling method and rename if necessary
                gfvs = list(cfg[section_name][label_cfg][gfm].keys())
                for gfv in gfvs:
                    if gfv in renames_exact:
                        new_name = std["Variables"]["rename_exact"][gfv]
                        cfg[section_name][label_cfg][gfm].rename(gfv, new_name)

                # now loop over the drivers and rename if necessary
                gfvs = list(cfg[section_name][label_cfg][gfm].keys())
                for gfv in gfvs:
                    drivers = cfg[section_name][label_cfg][gfm][gfv]["drivers"]
                    for re in renames_exact:
                        if re in drivers:
                            drivers = drivers.replace(re, std["Variables"]["rename_exact"][re])
                    cfg[section_name][label_cfg][gfm][gfv]["drivers"] = drivers

                # rename if first characters of variable name match pattern
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    # loop over the variables in the control file
                    gfvs = list(cfg[section_name][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        if (gfv.split("_")[0] == rp):
                            new_name = gfv.replace(gfv.split("_")[0], srp)
                            cfg[section_name][label_cfg][gfm].rename(gfv, new_name)
                # loop over the drivers
                gfvs = list(cfg[section_name][label_cfg][gfm].keys())
                for gfv in gfvs:
                    drivers = cfg[section_name][label_cfg][gfm][gfv]["drivers"]
                    for rp in renames_pattern:
                        if rp in drivers:
                            drivers = drivers.replace(rp, std["Variables"]["rename_pattern"][rp])
            else:
                continue

    return cfg

def update_cfg_variables_rename(cfg, std):
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
    renames_exact = list(std["Variables"]["rename_exact"].keys())
    # find out the section name
    for item in ["Variables", "Drivers", "Fluxes"]:
        if item in cfg:
            section_name = item
    # loop over the variables in the Variables section of the control file
    for label_cfg in list(cfg[section_name].keys()):
        # check and if required, rename the "name" key value
        if "name" in list(cfg[section_name][label_cfg].keys()):
            name = cfg[section_name][label_cfg]["name"]
            if name in renames_exact:
                new_name = std["rename_exact"][name]
                cfg[section_name][label_cfg]["name"] = new_name
        # check and if required, rename the variable itself
        if label_cfg in renames_exact:
            new_name = std["Variables"]["rename_exact"][label_cfg]
            renamed[label_cfg] = new_name
            cfg[section_name].rename(label_cfg, new_name)

    # rename pattern matches
    renames_pattern = list(std["Variables"]["rename_pattern"].keys())
    for rp in renames_pattern:
        #llen = len(rp)
        srp = std["Variables"]["rename_pattern"][rp]
        #klen = len(srp)
        # loop over the variables in the control file
        for label_cfg in list(cfg[section_name].keys()):
            # check and if required, rename the "name" key value
            if "name" in list(cfg[section_name][label_cfg].keys()):
                name = cfg[section_name][label_cfg]["name"]
                # 20210913 PRI - manual rebase
                #if ((name[:llen] == rp) and (name[:klen] != srp)):
                    #if len(name) > len(rp):
                        #if name[len(rp)] == "_":
                            #new_name = name.replace(name[:llen], srp)
                            #cfg[section_name][label_cfg]["name"] = new_name
                        #else:
                            ## different variable name, leave it alone
                            #pass
                    #else:
                        #new_name = name.replace(name[:llen], srp)
                        #cfg[section_name][label_cfg]["name"] = new_name
                if (name.split("_")[0] == rp):
                    new_name = name.replace(name.split("_")[0], srp)
                    cfg[section_name][label_cfg]["name"] = new_name
            # check and if required, rename the variable itself
            # 20210913 PRI - manual rebase
            #if ((label_cfg[:llen] == rp) and (label_cfg[:klen] != srp)):
                #if len(label_cfg) > len(rp):
                    #if label_cfg[len(rp)] == "_":
                        #new_name = label_cfg.replace(label_cfg[:llen], srp)
                        #renamed[label_cfg] = new_name
                        #cfg[section_name].rename(label_cfg, new_name)
                    #else:
                        ## different variable name, leave it alone
                        #pass
                #else:
                    #new_name = label_cfg.replace(label_cfg[:llen], srp)
                    #renamed[label_cfg] = new_name
                    #cfg[section_name].rename(label_cfg, new_name)
            if (label_cfg.split("_")[0] == rp):
                new_name = label_cfg.replace(label_cfg.split("_")[0], srp)
                renamed[label_cfg] = new_name
                cfg["Variables"].rename(label_cfg, new_name)

    # do any functions
    for label_cfg in list(cfg[section_name].keys()):
        if "Function" in cfg[section_name][label_cfg]:
            func_str = cfg[section_name][label_cfg]["Function"]["func"]
            for old_label in list(renamed.keys()):
                if old_label in func_str:
                    func_str = func_str.replace(old_label, renamed[old_label])
                    cfg[section_name][label_cfg]["Function"]["func"] = func_str

    # loop over the variables in the Plots section of the control file
    if "Plots" in list(cfg.keys()):
        plots = list(cfg["Plots"])
        for plot in plots:
            if "variables" in cfg["Plots"][plot]:
                cp = cfg["Plots"][plot]["variables"]
                vs = pfp_utils.string_to_list(cp)
                vs = [std["Variables"]["rename_exact"][v] if v in renames_exact else v for v in vs]
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    vs = [v.replace(v.split("_")[0], srp)
                          if (v.split("_")[0] == rp)
                          else v for v in vs]
                cfg["Plots"][plot]["variables"] = ",".join(vs)
            elif "type" in cfg["Plots"][plot]:
                if cfg["Plots"][plot]["type"] == "xy":
                    for axis in ["xseries", "yseries"]:
                        cp = cfg["Plots"][plot][axis]
                        vs = pfp_utils.string_to_list(cp)
                        vs = [std["Variables"]["rename_exact"][v] if v in renames_exact else v for v in vs]
                        for rp in renames_pattern:
                            srp = std["Variables"]["rename_pattern"][rp]
                            vs = [v.replace(v.split("_")[0], srp)
                                  if (v.split("_")[0] == rp)
                                  else v for v in vs]
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L3 control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_options(cfg, std)
        cfg = update_cfg_variables_deprecated(cfg, std)
        cfg = update_cfg_variables_rename(cfg, std)
        cfg = update_cfg_variable_attributes(cfg, std)
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

def update_cfg_options(cfg, std):
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
    if cfg["level"] == "L3":
        if "Options" in cfg:
            for item in cfg["Options"]:
                if item in ["ApplyFcStorage", "CcUnits", "FcUnits", "ReplaceFcStorage", "DisableFcWPL"]:
                    cfg["Options"].rename(item, item.replace("Fc", "Fco2"))
            old_units = list(std["Variables"]["units_map"].keys())
            for item in list(cfg["Options"].keys()):
                if cfg["Options"][item] in old_units:
                    cfg["Options"][item] = std["Variables"]["units_map"][cfg["Options"][item]]
    elif cfg["level"] == "L5":
        # rename exact variable name matches
        renames_exact = list(std["Variables"]["rename_exact"].keys())
        # rename pattern matches
        renames_pattern = list(std["Variables"]["rename_pattern"].keys())
        # loop over entries in the Options section, exact match
        options = cfg["Options"]
        for option in options:
            opt = cfg["Options"][option]
            if opt in renames_exact:
                cfg["Options"][option] = std["Variables"]["rename_exact"][opt]
        # loop over entries in the Options section, pattern match
        options = cfg["Options"]
        for rp in renames_pattern:
            llen = len(rp)
            srp = std["Variables"]["rename_pattern"][rp]
            klen = len(srp)
            for option in options:
                opt = str(cfg["Options"][option])
                if ((opt[:llen] == rp) and
                    (opt[:klen] != srp)):
                    new_opt = opt.replace(opt[:llen], srp)
                    cfg["Options"][option] = new_opt
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L4 control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_variables_deprecated(cfg, std)
        cfg = update_cfg_variables_rename(cfg, std)
        cfg = update_cfg_variable_attributes(cfg, std)
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
    # initialise the return logical
    ok = True
    try:
        cfg = update_cfg_syntax(cfg, std)
    except Exception:
        ok = False
        msg = " An error occurred while updating the L5 control file syntax"
    # clean up the variable names
    try:
        cfg = update_cfg_options(cfg, std)
        cfg = update_cfg_variables_deprecated(cfg, std)
        cfg = update_cfg_variables_rename(cfg, std)
        cfg = update_cfg_variable_attributes(cfg, std)
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
    # check to see if we can load the update_control_files.txt standard control file
    try:
        # get the base path of script or Pyinstaller application
        base_path = pfp_utils.get_base_path()
        stdname = os.path.join(base_path, "controlfiles", "standard", "update_control_files.txt")
        std = pfp_io.get_controlfilecontents(stdname, mode="quiet")
    except Exception:
        ok = False
        msg = " Unable to load standard control file " + stdname
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
    # clean up the variable names
    try:
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
    renames_exact = list(std["Variables"]["rename_exact"].keys())
    renames_pattern = list(std["Variables"]["rename_pattern"].keys())
    # loop over the variables in the EcosystemRespiration section of the control file
    for method in ["EcosystemRespiration"]:
        labels_cfg = list(cfg[method].keys())
        for label_cfg in labels_cfg:
            for qc in ["MergeSeries", "AverageSeries"]:
                if qc in cfg[method][label_cfg]:
                    source = cfg[method][label_cfg][qc]["source"]
                    vs = pfp_utils.string_to_list(source)
                    vs = [std["Variables"]["rename_exact"][v] if v in renames_exact else v for v in vs]
                    for rp in renames_pattern:
                        srp = std["Variables"]["rename_pattern"][rp]
                        vs = [v.replace(v.split("_")[0], srp)
                              if (v.split("_")[0] == rp)
                              else v for v in vs]
                    cfg[method][label_cfg][qc]["source"] = ",".join(vs)
                else:
                    continue
            for gfm in ["ERUsingSOLO", "ERUsingLloydTaylor", "ERUsingLasslop"]:
                if gfm in cfg[method][label_cfg]:
                    # loop over the variables for this gap filling method and rename if necessary
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        if gfv in renames_exact:
                            new_name = std["Variables"]["rename_exact"][gfv]
                            cfg[method][label_cfg][gfm].rename(gfv, new_name)

                    # now loop over the drivers and rename if necessary
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        item = cfg[method][label_cfg][gfm][gfv]["drivers"]
                        drivers = pfp_utils.string_to_list(item)
                        for driver in list(drivers):
                            for re in renames_exact:
                                if driver == re:
                                    new_name = std["Variables"]["rename_exact"][re]
                                    drivers[drivers.index(driver)] = new_name
                        cfg[method][label_cfg][gfm][gfv]["drivers"] = ",".join(drivers)

                    # rename if first characters of variable name match pattern
                    for rp in renames_pattern:
                        srp = std["Variables"]["rename_pattern"][rp]
                        # loop over the variables in the control file
                        gfvs = list(cfg[method][label_cfg][gfm].keys())
                        for gfv in gfvs:
                            if (gfv.split("_")[0] == rp):
                                new_name = gfv.replace(gfv.split("_")[0], srp)
                                cfg[method][label_cfg][gfm].rename(gfv, new_name)

                    # loop over the drivers
                    gfvs = list(cfg[method][label_cfg][gfm].keys())
                    for gfv in gfvs:
                        item = cfg[method][label_cfg][gfm][gfv]["drivers"]
                        drivers = pfp_utils.string_to_list(item)
                        for driver in list(drivers):
                            for rp in renames_pattern:
                                srp = std["Variables"]["rename_pattern"][rp]
                                if (driver.split("_")[0] == rp):
                                    new_name = opt.replace(driver.split("_")[0], srp)
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
                    new_name = std["Variables"]["rename_exact"][var]
                    cfg[method][label_cfg].rename(var, new_name)
                # pattern matches
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    if (var.split("_")[0] == rp):
                        new_name = var.replace(var.split("_")[0], srp)
                        cfg[method][label_cfg].rename(var, new_name)
            # second loop is over the values in case any of the keys were renamed
            variables = list(cfg[method][label_cfg].keys())
            for var in variables:
                opt = cfg[method][label_cfg][var]
                # exact matches
                if opt in renames_exact:
                    new_name = std["Variables"]["rename_exact"][opt]
                    cfg[method][label_cfg][var] = new_name
                # pattern matches
                for rp in renames_pattern:
                    srp = std["Variables"]["rename_pattern"][rp]
                    if (opt.split("_")[0] == rp):
                        new_name = opt.replace(opt.split("_")[0], srp)
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
    renames_exact = list(std["Variables"]["rename_exact"].keys())
    renames_pattern = list(std["Variables"]["rename_pattern"].keys())
    # loop over the variables in the EcosystemRespiration, NetEcosystemExchange and
    # GrossPrimaryProductivity section of the control file
    for method in ["EcosystemRespiration", "NetEcosystemExchange", "GrossPrimaryProductivity"]:
        # rename exact matches
        for label_cfg in list(cfg[method].keys()):
            if label_cfg in renames_exact:
                new_name = std["Variables"]["rename_exact"][label_cfg]
                cfg[method].rename(label_cfg, new_name)
        # rename pattern matches
        for rp in renames_pattern:
            srp = std["Variables"]["rename_pattern"][rp]
            # loop over the variables in the control file
            for label_cfg in list(cfg[method].keys()):
                if (label_cfg.split("_")[0] == rp):
                    new_name = label_cfg.replace(label_cfg.split("_")[0], srp)
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
    if pfp_utils.is_number(v):
        pass
    elif ("[" in v) and ("]" in v) and ("*" in v):
        # old style of [value]*12
        v = v[v.index("[")+1:v.index("]")]
    elif ("[" in v) and ("]" in v) and ("*" not in v):
        # old style of [1,2,3,4,5,6,7,8,9,10,11,12]
        v = v.replace("[", "").replace("]", "")
    else:
        pass
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
