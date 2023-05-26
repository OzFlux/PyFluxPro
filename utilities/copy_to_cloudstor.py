# standard modules
import logging
import os
import subprocess
import sys
import time
# 3rd party modules
import xlrd
# PFP modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("utilities")-1])
import scripts.pfp_io as pfp_io
import scripts.pfp_utils as pfp_utils

logger = logging.getLogger("copy_to_cloudstor")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
console = logging.StreamHandler()
console.setFormatter(formatter)
console.setLevel(logging.DEBUG)
logger.addHandler(console)

def check_csep_folder(csep_path, create=True, dry_run=True):
    ok = True
    rclone_cmd = ["rclone", "lsd", csep_path]
    if dry_run:
        msg = " ".join(rclone_cmd)
        logger.info(msg)
        return_code = 0
    else:
        return_code = subprocess.call(rclone_cmd,
                                      stdout=subprocess.DEVNULL,
                                      stderr=subprocess.STDOUT)
    if return_code == 0:
        pass
    elif return_code == 3 and create:
        rclone_cmd = ["rclone", "mkdir", csep_path]
        if dry_run:
            msg = " ".join(rclone_cmd)
            logger.info(msg)
            return_code = 0
        else:
            return_code = subprocess.call(rclone_cmd,
                                          stdout=subprocess.DEVNULL,
                                          stderr=subprocess.STDOUT)
        if return_code != 0:
            msg = "An error occurred creating directory " + csep_path
            logger.info(msg)
            ok = False
    else:
        msg = "Unexpected return code " + str(return_code)
        logger.info(msg)
        ok = False
    return ok
def copy_to_remote(source_name, csep_uri, dry_run=True):
    ok = True
    start = time.time()
    rclone_cmd = ["rclone", "copyto", source_name, csep_uri]
    if dry_run:
        msg = " ".join(rclone_cmd)
        logger.info(msg)
        return_code = 0
    else:
        return_code = subprocess.call(rclone_cmd,
                                      stdout=subprocess.DEVNULL,
                                      stderr=subprocess.STDOUT)
    if return_code != 0:
        msg = "An error occurred copying " + source_name
        logger.error(msg)
        ok = False
    else:
        end = time.time()
        elapsed_time = end - start
        MBps = os.stat(source_name).st_size/float(1024*1024*elapsed_time)
        msg = " Copied " + source_name + " (" + str(round(MBps, 3)) + " MB/s)"
        logger.info(msg)
    return ok
def get_destination_name(source_name):
    ok = True
    destination_name = ""
    if not os.path.isfile(source_name):
        msg = " File " + source_name + " not found"
        logger.error(msg)
        ok = False
    else:
        ds = pfp_io.NetCDFRead(source_name)
        ldt = pfp_utils.GetVariable(ds, "DateTime")
        start_datetime = ldt["Data"][0]
        end_datetime = ldt["Data"][-1]
        file_datetime = "_" + start_datetime.strftime("%Y%m%d") + "_"
        file_datetime += end_datetime.strftime("%Y%m%d")
        destination_name = source_name.replace(".nc", file_datetime+".nc")
    return destination_name, ok

logger.info("Reading site_master.xls")
site_master_name = "/mnt/OzFlux/Sites/site_master.xls"
xlwb = xlrd.open_workbook(site_master_name)
processing_sheet = xlwb.sheet_by_name("Processing")
processing_info = {}
labels = [h for h in processing_sheet.row_values(0)]
for n, label in enumerate(labels):
    processing_info[label] = processing_sheet.col_values(n)[1:]
#bp_site_names = sorted(processing_info["BP name"])
bp_site_names = ["Whroo"]
#msg = "Processing " + ",".join(bp_site_names)
#logger.info(msg)

bp_base = "/mnt/OzFlux/Sites"
csep_base = "TERNcloudstor:Shared/EcoSystemProcesses/Sites"
csep_version = "2023_v1"
levels = ["L3", "L4", "L5", "L6"]
method = "default"
create = True
dry_run = False

for bp_site_name in bp_site_names:
    csep_path = csep_base
    idx = processing_info["BP name"].index(bp_site_name)
    csep_site_name = processing_info["CSEP name"][idx]
    csep_folders = [csep_site_name, "Data", "Flux", "Processed", csep_version]
    for csep_folder in csep_folders:
        csep_path = os.path.join(csep_path, csep_folder)
        if not check_csep_folder(csep_path, create=create, dry_run=dry_run):
            continue
    csep_version_path = csep_path
    for level in levels:
        csep_level_path = os.path.join(csep_version_path, level)
        if not check_csep_folder(csep_level_path, create=create, dry_run=dry_run):
            continue
        csep_method_path = os.path.join(csep_level_path, method)
        if not check_csep_folder(csep_method_path, create=create, dry_run=dry_run):
            continue
        source_path = os.path.join(bp_base, bp_site_name, "Data", "Processed")
        os.chdir(source_path)
        source_name = bp_site_name + "_" + level + ".nc"
        destination_name, ok = get_destination_name(source_name)
        if not ok:
            continue
        csep_path = os.path.join(csep_base, csep_site_name, "Data", "Flux", "Processed")
        csep_uri = os.path.join(csep_method_path, destination_name)
        if not copy_to_remote(source_name, csep_uri, dry_run=dry_run):
            continue
        if level == "L6":
            for suffix in ["_Summary", "_Annual", "_Cumulative", "_Daily", "_Monthly"]:
                summary_name = source_name.replace(".nc", suffix+".nc")
                if not os.path.isfile(summary_name):
                    msg = " File " + summary_name + " not found"
                    logger.error(msg)
                    continue
                summary_destination_name = destination_name.replace(".nc", suffix+".nc")
                csep_uri = os.path.join(csep_method_path, summary_destination_name)
                if not copy_to_remote(summary_name, csep_uri, dry_run=dry_run):
                    continue
msg = "Finished"
logger.info(msg)