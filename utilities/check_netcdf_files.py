# standard modules
import datetime
import logging
import os
import subprocess
import sys
# 3rd party modules
from configobj import ConfigObj
import numpy
import timezonefinder
# PFP modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("utilities")-1])
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_utils as pfp_utils

now = datetime.datetime.now()
log_file_name = "check_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_path = "logfiles"
if not os.path.isdir(log_file_path):
    os.makedirs(log_file_path)
log_file_name = os.path.join("logfiles", log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_name, to_screen=True)
cfchecker_file_name = log_file_name.replace(".log", "_cfchecker.all")
cfchecker_file = open(cfchecker_file_name, "w")
error_file_name = cfchecker_file_name.replace(".all", ".errors")
warning_file_name = cfchecker_file_name.replace(".all", ".warnings")

def check_global_attributes(ds, chk, messages):
    # check the 'Global' section exists
    if not hasattr(ds, "globalattributes"):
        msg = "Global: no global attributes section"
        messages["ERROR"].append(msg)
        return
    if len(ds.globalattributes.keys()) == 0:
        msg = "Global: no global attributes defined"
        messages["ERROR"].append(msg)
        return
    # check the required global attributes exist
    check_global_required(ds, chk, messages)
    # check the forced global attributes
    check_global_forced(ds, chk, messages)
    # check the recommended global attributes
    check_global_recommended(ds, chk, messages)
    # check the deprecated global attributes
    check_global_deprecated(ds, chk, messages)
    return
def check_global_required(ds, chk, messages):
    # check the global attributes
    required = chk["Global"]["Required"]
    ds_global = ds.globalattributes
    # check the required global attributes are present
    for item in required:
        if item not in ds_global:
            msg = "Global: " + item + " not in section (required)"
            messages["ERROR"].append(msg)
    # check time step is present and makes sense
    if "time_step" in ds_global:
        try:
            ts = int(float(ds_global["time_step"]))
        except ValueError:
            msg = "Global: 'time_step' is not a number"
            messages["ERROR"].append(msg)
        if ts not in [15, 20, 30, 60]:
            msg = "Global : 'time_step' must be 15, 20, 30 or 60"
            messages["ERROR"].append(msg)
    # check latitude is present and makes sense
    if "latitude" in ds_global:
        try:
            lat = float(ds_global["latitude"])
            if lat < -90.0 or lat > 90.0:
                msg = "Global: 'latitude' must be between -90 and 90"
                messages["ERROR"].append(msg)
        except ValueError:
            msg = "Global: 'latitude' is not a number"
            messages["ERROR"].append(msg)
    # check longitude is present and makes sense
    if "longitude" in ds_global:
        try:
            lon = float(ds_global["longitude"])
            if lon < -180.0 or lat > 180.0:
                msg = "Global: 'longitude' must be between -180 and 180"
                messages["ERROR"].append(msg)
        except ValueError:
            msg = "Global: 'longitude' is not a number"
            messages["ERROR"].append(msg)
    return
def check_global_forced(ds, chk, messages):
    forced = chk["Global"]["Forced"]
    ds_global = ds.globalattributes
    # force global attributes that have defined values
    for item in forced:
        if item not in ds_global:
            msg = "Global: " + item + " not found in global attributes"
            messages["ERROR"].append(msg)
            continue
        else:
            if ds_global[item] != forced[item]:
                msg = "Global: " + item + " (" + ds_global[item] + "), should be " + forced[item]
                messages["ERROR"].append(msg)
    # and do the time zone
    lon = float(ds_global["longitude"])
    lat = float(ds_global["latitude"])
    tf = timezonefinder.TimezoneFinder()
    time_zone = tf.timezone_at(lng=lon, lat=lat)
    if "time_zone" in ds_global:
        if ds_global["time_zone"] != time_zone:
            msg = "Global: existing time zone (" + ds_global["time_zone"] + ") doesn't match " + time_zone
            messages["ERROR"].append(msg)
        else:
            pass
    else:
        ds_global["time_zone"] = time_zone
    return
def check_global_recommended(ds, chk, messages):
    recommended = chk["Global"]["Recommended"]
    ds_global = ds.globalattributes
    # check recommended global attributes
    for item in recommended:
        if item not in ds_global:
            msg = "Global: recommended global attribute " + item + " not found"
            messages["ERROR"].append(msg)
    return
def check_global_deprecated(ds, chk, messages):
    deprecated = chk["Global"]["Deprecated"]
    ds_global = ds.globalattributes
    for item in deprecated:
        if item in ds_global:
            msg = "Global: " + item + " has been deprecated"
            if len(deprecated[item]) != 0:
                msg += ", replace with " + deprecated[item]
            messages["ERROR"].append(msg)
    return
def check_variables_sections(ds, chk, ds_label, chk_label, messages):
    var_keys = list(ds.series[ds_label].keys())
    # check we have an 'Attr' section, remove variable if absent
    if ("Attr" in var_keys):
        var_attrs = list(ds.series[ds_label]["Attr"].keys())
        if "long_name" in var_attrs:
            pass
        else:
            msg = "Variables: " + ds_label + "; no long_name variable attribute"
            messages["ERROR"].append(msg)
        # check statistic_type
        if "statistic_type" in var_attrs:
            check_variable_statistic_type(ds, chk, ds_label, chk_label, messages)
            # check units
            if "units" in var_attrs:
                check_variable_units(ds, chk, ds_label, chk_label, messages)
                check_variable_attributes_consistent(ds, chk, ds_label, chk_label, messages)
                if "standard_name" in var_attrs:
                    check_variable_standard_name(ds, chk, ds_label, chk_label, messages)
                else:
                    pass
            else:
                msg = "Variables: " + ds_label + "; no units variable attribute"
                messages["ERROR"].append(msg)
        else:
            msg = "Variables: " + ds_label + "; no statistic_type variable attribute"
            messages["ERROR"].append(msg)
        # check height given for CO2 value
        if "CO2" in ds_label:
            check_variable_height(ds, ds_label, messages)
        else:
            pass
    else:
        msg = "Variables: " + ds_label + "; 'Attr' section missing"
        messages["ERROR"].append(msg)
    return
def check_variable_statistic_type(ds, chk, ds_label, chk_label, messages):
    ds_attr = ds.series[ds_label]["Attr"]
    chk_var = chk["Variables"][chk_label]
    statistic_types = sorted(list(chk_var.keys()))
    ds_stat_type = ds_attr["statistic_type"]
    if ds_stat_type not in statistic_types:
        msg = "Variables: " + ds_label + "; unrecognised statistic_type (" + ds_stat_type + ")"
        messages["ERROR"].append(msg)
    return
def check_variable_units(ds, chk, ds_label, chk_label, messages):
    ds_attr = ds.series[ds_label]["Attr"]
    ds_stat_type = ds_attr["statistic_type"]
    if ds_stat_type in chk["Variables"][chk_label]:
        chk_stat_type = chk["Variables"][chk_label][ds_stat_type]
        chk_units = sorted(list(chk_stat_type.keys()))
        ds_units = ds_attr["units"]
        if ds_units not in chk_units:
            msg = "Variables: " + ds_label + "; unrecognised units (" + ds_units + ")"
            messages["ERROR"].append(msg)
    else:
        msg = "Variables: " + ds_label + "; unrecognised statistic type (" + ds_stat_type + ")"
        messages["ERROR"].append(msg)
    return
def check_variable_attributes_consistent(ds, chk, ds_label, chk_label, messages):
    ds_attr = ds.series[ds_label]["Attr"]
    ds_units = ds_attr["units"]
    ds_stat_type = ds_attr["statistic_type"]
    if ds_stat_type in chk["Variables"][chk_label]:
        chk_stat_type = chk["Variables"][chk_label][ds_stat_type]
        if ds_units in chk_stat_type:
            for item in chk_stat_type[ds_units]:
                if (item not in ds_attr):
                    # attribute not found
                    msg = "Variables: " + ds_label + "; attribute (" + item + ") not found"
                    messages["ERROR"].append(msg)
                elif (ds_attr[item] != chk_stat_type[ds_units][item]):
                    # attribute not the same as the chk attribute
                    msg = "Variables: " + ds_label + "; invalid " + item + " (" + ds_attr[item] + ")"
                    messages["ERROR"].append(msg)
                else:
                    pass
        else:
            msg = "Variables: " + ds_label + "; unrecognised units (" + ds_units + ")"
            if msg not in messages["ERROR"]:
                messages["ERROR"].append(msg)
    else:
        msg = "Variables: " + ds_label + "; unrecognised statistic type (" + ds_stat_type + ")"
        if msg not in messages["ERROR"]:
            messages["ERROR"].append(msg)
    return
def check_variable_standard_name(ds, chk, ds_label, chk_label, messages):
    ds_attr = ds.series[ds_label]["Attr"]
    ds_units = ds_attr["units"]
    ds_stat_type = ds_attr["statistic_type"]
    if ds_stat_type in chk["Variables"][chk_label]:
        chk_stat_type = chk["Variables"][chk_label][ds_stat_type]
        if ds_units in chk_stat_type:
            if (("standard_name" in ds_attr) and
                ("standard_name" not in chk_stat_type[ds_units])):
                msg = "Variables: " + ds_label + "; standard_name not allowed"
                messages["ERROR"].append(msg)
            else:
                pass
        else:
            msg = "Variables: " + ds_label + "; unrecognised units (" + ds_units + ")"
            if msg not in messages["ERROR"]:
                messages["ERROR"].append(msg)
    else:
        msg = "Variables: " + ds_label + "; unrecognised statistic type (" + ds_stat_type + ")"
        if msg not in messages["ERROR"]:
            messages["ERROR"].append(msg)
    return
def check_variable_height(ds, ds_label, messages):
    ds_attr = ds.series[ds_label]["Attr"]
    if "height" in ds_attr:
        height = pfp_utils.strip_non_numeric(ds_attr["height"])
        if pfp_utils.is_number(height):
            pass
        else:
            msg = "Variables: " + ds_label + "; 'height' attribute not a number"
            messages["ERROR"].append(msg)
    else:
        msg = "Variables: " + ds_label + "; no 'height' attribute found"
        messages["ERROR"].append(msg)
    return
def CheckCFCompliance(nc_file_uri, messages):
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
        # with a little work, cfchecker could be used directly in PFP
        # path to Python running this show
        env_bin = os.path.join(sys.exec_prefix, "bin")
        # path to cfchecks in that Python install
        env_cfchecks = os.path.join(env_bin, "cfchecks")
        # command to run with path to cfchecks
        cmd = [env_cfchecks, "-v 1.8", nc_file_uri]
        with open(cf_file_uri, "w") as f:
            # run cfchecks as a subprocess at this stage
            subprocess.run(cmd, stdout=f)
        # parse the output of cfchecker and write warnings and errors to the message dictionary
        # read the cfchecker log file
        with open(cf_file_uri) as f:
            for line in f:
                if "Checking variable:" in line:
                    # line containing the variable name
                    parts = line.split(" ")
                    variable_name = parts[2].rstrip()
                    continue
                if (("ERROR:" in line) and ('variable_name' in locals())):
                    # line containing error messages
                    msg = "CF: " + variable_name + "; " + line.rstrip()
                    messages["ERROR"].append(msg)
                    continue
                if (("WARN:" in line) and ('variable_name' in locals())):
                    # line containing warning messages
                    msg = "CF: " + variable_name + "; " + line.rstrip()
                    messages["WARNING"].append(msg)
                    continue
    except Exception:
        msg = "Error during CF compliance check of " + os.path.split(nc_file_uri)[1]
        logger.error(msg)
    return

logger = logging.getLogger("pfp_log")
base_path = "/mnt/OzFlux/Sites"
site_names = ["Whroo"]
#site_names = ["AdelaideRiver", "AliceSpringsMulga", "Boyagin", "Calperum", "CapeTribulation", "Collie",
              #"CowBay", "CumberlandPlain", "DalyPasture", "DalyUncleared", "DryRiver", "Emerald",
              #"FoggDam", "Gingin", "GreatWesternWoodlands", "HowardSprings", "Litchfield", "Longreach",
              #"Loxton", "Otway", "RedDirtMelonFarm", "Ridgefield", "RiggsCreek", "RobsonCreek", "Samford",
              #"SilverPlains", "SturtPlains", "TiTreeEast", "Tumbarumba", "WallabyCreek", "Warra", "Whroo",
              #"WombatStateForest", "Yanco"]

chk_name = "/home/pisaac/PyFluxPro/controlfiles/standard/check_l1_controlfile.txt"
chk = ConfigObj(chk_name, indent_type="    ", list_values=False, write_empty_values=True)
chk_labels = sorted(list(chk["Variables"].keys()))

for site_name in site_names:
    msg = " Processing site " + site_name
    logger.info(msg)
    site_portal_path = os.path.join(base_path, site_name, "Data", "Processed")
    file_names = sorted([f.name for f in os.scandir(site_portal_path) if f.is_file()])
    file_names = [f for f in file_names if (("L3" in f) and (".nc" in f))]
    msg = "  Files: " + ",".join(file_names)
    logger.info(msg)
    for file_name in file_names:
        nc_file_uri = os.path.join(site_portal_path, file_name)
        # read the netCDF file, don't correct the time step or update the file
        ds = pfp_io.NetCDFRead(nc_file_uri, checktimestep=False, update=False)
        ldt = pfp_utils.GetVariable(ds, "DateTime")
        ts = int(float(ds.globalattributes["time_step"]))
        dt = pfp_utils.get_timestep(ds)
        index = numpy.where(dt != ts*60)[0]
        if len(index) != 0:
            msg = str(len(index)) + " problems found with the time stamp"
            logger.warning(msg)
            msg = "The first 10 are:"
            logger.warning(msg)
            for i in range(min([10, len(index)])):
                msg = "  " + str(ldt["Data"][i-1]) + str(ldt["Data"][i]) + str(ldt["Data"][i+1])
                logger.warning(msg)
        # read the netCDF file again, this time correct the time step but don't update the file
        ds = pfp_io.NetCDFRead(nc_file_uri, update=False)
        pyfluxpro_version = ""
        for item in ["QC_version", "pyfluxpro_version"]:
            if item in ds.globalattributes:
                pyfluxpro_version = ds.globalattributes[item]
                break
        for item in ["nc_rundatetime", "date_created"]:
            if item in ds.globalattributes:
                date_created = ds.globalattributes[item]
                break
        msg = " Processed using " + pyfluxpro_version + " on " + date_created
        logger.info(msg)
        ds_labels = sorted(list(ds.series.keys()))
        # initialise the messages dictionary
        messages = {"ERROR":[], "WARNING": [], "INFO": []}
        # check the global attributes section
        check_global_attributes(ds, chk, messages)
        # check variables whose name exactly matches an entry in the settings/l1.txt control file
        done = []
        label_matches = [l for l in ds_labels if l in chk_labels]
        for ds_label in label_matches:
            chk_label = ds_label
            # check variable 'Attr' section
            check_variables_sections(ds, chk, ds_label, chk_label, messages)
            # append this variable name to the done list
            done.append(ds_label)
        # check variables where the first characters of the name match an entry in settings/l1.txt
        ds_labels = sorted(list(ds.series.keys()))
        for chk_label in chk_labels:
            lcl = len(chk_label)
            label_matches = [l for l in ds_labels if l.split("_")[0] == chk_label and l not in done]
            for ds_label in label_matches:
                # check variable 'Attr' section
                check_variables_sections(ds, chk, ds_label, chk_label, messages)
                # append this variable name to the done list
                done.append(ds_label)
        
        CheckCFCompliance(nc_file_uri, messages)
        
        #print(os.path.basename(nc_file_uri))
        for e in messages:
            #print("  ", e)
            for m in messages[e]:
                if "ERROR" in e:
                    logger.error("    " + m)
                elif "WARNING" in e:
                    logger.warning("    " + m)
                elif "INFO" in e:
                    logger.info("    " + m)
    logger.info("")