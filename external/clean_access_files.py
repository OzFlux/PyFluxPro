# standard modules
import datetime
import faulthandler
#import logging
import os
import sys
# 3rd party modules
import numpy
import xlrd
# PFP modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("external")-1])
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_utils as pfp_utils

faulthandler.enable()

## log to screen
#logToScreen = logging.getLogger("logToScreen")
#logToScreen.setLevel(logging.DEBUG)
#sh = logging.StreamHandler()
#sh.setLevel(logging.DEBUG)
#formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
#sh.setFormatter(formatter)
#logToScreen.addHandler(sh)

now = datetime.datetime.now()
log_file_name = "cleanup_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_path = "logfiles"
if not os.path.isdir(log_file_path):
    os.mkdir(log_file_path)
log_file_uri = os.path.join(log_file_path, log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_uri, to_screen=True)

def read_site_master(site_master_file_path, sheet_name):
    """
    """
    xl_book = xlrd.open_workbook(site_master_file_path)
    xl_sheet = xl_book.sheet_by_name(sheet_name)
    last_row = int(xl_sheet.nrows)
    # find the header and first data rows
    for i in range(last_row):
        if xl_sheet.cell(i,0).value == "Site":
            header_row = i
            first_data_row = header_row + 1
            break
    # read the header row
    header_row_values = xl_sheet.row_values(header_row)
    # read the site data from the master Excel spreadsheet
    site_info = {}
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = {}
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value
    return site_info

cfg_file_path = "ACCESS.txt"
msg = " Loading the control file"
#logToScreen.info(msg)
logger.info(msg)
cfg = pfp_io.get_controlfilecontents(cfg_file_path)
#cfg_labels = [l for l in list(cfg["Variables"].keys()) if "nc" in list(cfg["Variables"][l].keys())]
# read the site master workbook
site_info = read_site_master(cfg["Files"]["site_master_file_path"], "ACCESS")
sites = list(site_info.keys())
#sites = ["Ridgefield"]
n_sites = len(sites)
base_path = cfg["Files"]["existing_access_base_path"]

remove_pattern = ["Ah", "q", "Fn_sw", "Fn_lw"]
rename_pattern = {"u_": "U_", "v_": "V_"}

# loop over sites
for site in sites:
    #logToScreen.info("Processing " + site)
    logger.info(msg)
    access_name = site + "_ACCESS.nc"
    access_uri = os.path.join(base_path, site, "Data", "ACCESS", "downloaded", access_name)
    if not os.path.isfile(access_uri):
        msg = "ACCESS file " + access_name + " not found"
        #logToScreen.info(msg)
        logger.info(msg)
        continue
    ds = pfp_io.NetCDFRead(access_uri)
    # remove unwanted variables
    labels = list(ds.series.keys())
    for label in labels:
        for rp in remove_pattern:
            if label[:len(rp)] == rp:
                ds.series.pop(label)

    # rename variables
    labels = list(ds.series.keys())
    for rp_label in list(rename_pattern.keys()):
        access_labels = [l for l in labels if l[:len(rp_label)] == rp_label]
        for access_label in access_labels:
            new_label = access_label.replace(rp_label, rename_pattern[rp_label])
            ds.series[new_label] = ds.series.pop(access_label)

    # remove unwanted variable attributes
    labels = list(ds.series.keys())
    cfg_labels = [l for l in list(cfg["Variables"].keys())]
    for cfg_label in cfg_labels:
        access_labels = sorted([l for l in labels if l[:len(cfg_label)] == cfg_label])
        for access_label in access_labels:
            var = pfp_utils.GetVariable(ds, access_label)
            vattrs = list(var["Attr"].keys())
            for vattr in vattrs:
                if vattr not in list(cfg["Variables"][cfg_label]["Attr"].keys()):
                    var["Attr"].pop(vattr)
            pfp_utils.CreateVariable(ds, var)

    # force the remaining attributes to values in control file
    labels = list(ds.series.keys())
    cfg_labels = [l for l in list(cfg["Variables"].keys())]
    for cfg_label in cfg_labels:
        access_labels = [l for l in labels if l[:len(cfg_label)] == cfg_label]
        for access_label in access_labels:
            var = pfp_utils.GetVariable(ds, access_label)
            for vattr in list(cfg["Variables"][cfg_label]["Attr"]):
                var["Attr"][vattr] = cfg["Variables"][cfg_label]["Attr"][vattr]
            pfp_utils.CreateVariable(ds, var)

    # force missing data to non-zero QC flag
    labels = list(ds.series.keys())
    for label in labels:
        var = pfp_utils.GetVariable(ds, label)
        condition = numpy.ma.getmaskarray(var["Data"]) & (numpy.mod(var["Flag"],10) == 0)
        idx = numpy.ma.where(condition == True)[0]
        if len(idx) != 0:
            #msg = " "+label+": "+str(len(idx))+" missing values with flag = 0 (forced to 8)"
            #logger.warning(msg)
            var["Flag"][idx] = numpy.int32(8)
            pfp_utils.CreateVariable(ds, var)

    # rename the new ACCESS file
    #uid = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #access_rename_uri = access_uri.replace(".nc", "_"+uid+".nc")
    #os.rename(access_uri, access_rename_uri)
    # write out the new ACCESS file
    msg = "Writing " + access_name
    #logToScreen.info(msg)
    logger.info(msg)
    #access_path = os.path.join(base_path, site, "Data", "ACCESS", "cleaned")
    access_path = os.path.join(base_path, site, "Data", "ACCESS")
    if not os.path.isdir(access_path):
        os.mkdir(access_path)
    access_uri = os.path.join(access_path, access_name)
    pfp_io.NetCDFWrite(access_uri, ds)
#logToScreen("Finished")
logger.info("Finished")