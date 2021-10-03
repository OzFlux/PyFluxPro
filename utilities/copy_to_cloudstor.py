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

bp_base = "/mnt/OzFlux/Sites"
csep_base = "TERNcloudstor:Shared/EcoSystemProcesses/Sites"

logger.info("Reading site_master.xls")
site_master_name = "/mnt/OzFlux/Sites/site_master_pri.xls"
xlwb = xlrd.open_workbook(site_master_name)
processing_sheet = xlwb.sheet_by_name("Processing")
processing_info = {}
labels = [h for h in processing_sheet.row_values(0)]
for n, label in enumerate(labels):
    processing_info[label] = processing_sheet.col_values(n)[1:]

bp_site_names = sorted(processing_info["BP name"])
msg = "Processing " + ",".join(bp_site_names)
logger.info(msg)
#bp_site_names = ["AdelaideRiver"]
for bp_site_name in bp_site_names:
    idx = processing_info["BP name"].index(bp_site_name)
    csep_site_name = processing_info["CSEP name"][idx]
    source_path = os.path.join(bp_base, bp_site_name, "Data", "Processed")
    os.chdir(source_path)
    msg = "In " + os.getcwd()
    logger.info(msg)
    for level in ["L3", "L4", "L5", "L6"]:
        source_name = bp_site_name + "_" + level + ".nc"
        if not os.path.isfile(source_name):
            msg = " File " + source_name + " not found"
            logger.error(msg)
            continue
        ds = pfp_io.NetCDFRead(source_name)
        ldt = pfp_utils.GetVariable(ds, "DateTime")
        start_datetime = ldt["Data"][0]
        end_datetime = ldt["Data"][-1]
        file_datetime = "_" + start_datetime.strftime("%Y%m%d") + "_"
        file_datetime += end_datetime.strftime("%Y%m%d")
        destination_name = source_name.replace(".nc", file_datetime+".nc")
        csep_path = os.path.join(csep_base, csep_site_name, "Data", "Flux", "Processed")
        csep_uri = os.path.join(csep_path, level, "default", destination_name)
        msg = " Copying " + source_name + " to " + destination_name
        logger.info(msg)
        start = time.time()
        rclone_cmd = ["rclone", "copyto", source_name, csep_uri]
        subprocess.call(rclone_cmd)
        #logger.info(" ".join(rclone_cmd))
        #time.sleep(1)
        end = time.time()
        elapsed_time = end - start
        MBps = os.stat(source_name).st_size/float(1024*1024*elapsed_time)
        msg = " Copied " + source_name + " (" + str(round(MBps, 3)) + " MB/s)"
        logger.info(msg)
        if level == "L6":
            for suffix in ["_Summary", "_Annual", "_Cumulative", "_Daily", "_Monthly"]:
                summary_name = source_name.replace(".nc", suffix+".nc")
                if not os.path.isfile(summary_name):
                    msg = " File " + summary_name + " not found"
                    logger.error(msg)
                    continue
                summary_destination_name = destination_name.replace(".nc", suffix+".nc")
                csep_uri = os.path.join(csep_path, level, "default", summary_destination_name)
                msg = " Copying " + summary_name + " to " + summary_destination_name
                logger.info(msg)
                start = time.time()
                rclone_cmd = ["rclone", "copyto", summary_name, csep_uri]
                subprocess.call(rclone_cmd)
                #logger.info(" ".join(rclone_cmd))
                #time.sleep(1)
                end = time.time()
                elapsed_time = end - start
                MBps = os.stat(summary_name).st_size/float(1024*1024*elapsed_time)
                msg = " Copied " + summary_name + " (" + str(round(MBps, 3)) + " MB/s)"
                logger.info(msg)
msg = "Finished"
logger.info(msg)