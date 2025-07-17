import datetime
import logging
from multiprocessing import Pool
import os
import sys
import warnings
# 3rd party modules
from configobj import ConfigObj
# PFP modules
sys.path.append('scripts')
from scripts import pfp_batch
from scripts import pfp_io
from scripts import pfp_log

warnings.filterwarnings("ignore", category=Warning)

logger = logging.getLogger("pfp_log")

def do_sites_batch(main_ui):
    """
    Purpose:
     Function to farm out batch processing by site across multiple CPUs.
     This function controls the multiprocessing, do_sites_batch_dispatcher
     does the work.

     Processing a site is independent of processing any other site.  This means
     that we can use multiprocessing to processes as many sites as possible
     in paralell because they will not interfere with each other.  This can't be
     done for batch processing by levels since the levels must be performed in
     the correct order.
    Usage:
    Side effects:
    Author: PRI
    Date: April 2022
    """
    # list to hold arguments for do_sites_batch_dispatcher
    cfs = []
    # loop over the sites
    for site in list(main_ui.cfg["Sites"].keys()):
        # construct the argument for do_sites_batch_dispatcher for this site
        obj = pfp_batch.Bunch(stop_flag=False, cfg=main_ui.cfg, mode="batch", site=site)
        # add it to the list
        cfs.append(obj)
    #for item in cfs:
        #do_sites_batch_dispatcher(item)
    # get the number of CPUs to use, this is the minimum of the number available
    # minus 1 (keep 1 CPU for foreground tasks) and 10
    number_cpus = min([os.cpu_count()-1, 10])
    # let Pool do its thing
    with Pool(number_cpus) as pool:
        pool.map(do_sites_batch_dispatcher, cfs)
    return

def do_sites_batch_dispatcher(item):
    """
    Purpose:
     This function loops over the control files for a given site listed in the
     batch control file, opens the control file to get the processing level
     and then calls the appropriate processing routine.
    Usage:
    Side effects:
    Author: PRI
    Date: April 2022
    """
    # local pointers
    # get the control files for this site
    cfg_site = item.cfg["Sites"][item.site]
    # loop over the control files
    for n in sorted(list(cfg_site.keys()), key=int):
        # get the control file name
        cfg_name = cfg_site[n]
        # get the control file contents
        cfg = pfp_io.get_controlfilecontents(cfg_name)
        # get the processing level
        level = str(cfg["level"])
        # get the batch routine argument, this is a dictionary with a single
        # entry containing the name of the control file to process
        cf_level = {n: item.cfg["Sites"][item.site][n]}
        # call the batch routine based on the processing level
        if level.lower() == "l1":
            # L1 processing
            pfp_batch.do_L1_batch(item, cf_level)
        elif level.lower() == "l2":
            # L2 processing
            pfp_batch.do_L2_batch(item, cf_level)
        elif level.lower() == "l3":
            # L3 processing
            pfp_batch.do_L3_batch(item, cf_level)
        elif level.lower() == "concatenate":
            # concatenate netCDF files
            pfp_batch.do_concatenate_batch(item, cf_level)
        elif level.lower() == "climatology":
            # climatology
            pfp_batch.do_climatology_batch(item, cf_level)
        elif level.lower() == "cpd_barr":
            # ustar threshold from change point detection
            pfp_batch.do_cpd_barr_batch(item, cf_level)
        elif level.lower() == "cpd_mchugh":
            # ustar threshold from change point detection
            pfp_batch.do_cpd_mchugh_batch(item, cf_level)
        elif level.lower() == "cpd_mcnew":
            # ustar threshold from change point detection
            pfp_batch.do_cpd_mcnew_batch(item, cf_level)
        elif level.lower() == "mpt":
            # ustar threshold from change point detection
            pfp_batch.do_mpt_batch(item, cf_level)
        elif level.lower() == "l4":
            # L4 processing
            pfp_batch.do_L4_batch(item, cf_level)
        elif level.lower() == "l5":
            # L5 processing
            pfp_batch.do_L5_batch(item, cf_level)
        elif level.lower() == "l6":
            # L6 processing
            pfp_batch.do_L6_batch(item, cf_level)
        else:
            msg = " Unrecognised batch processing level " + str(level)
            logger.error(msg)

if (__name__ == '__main__'):

    # get the logger
    now = datetime.datetime.now()
    log_file_name = "pfp_" + now.strftime("%Y%m%d%H%M") + ".log"
    log_file_name = os.path.join("logfiles", log_file_name)
    logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_name,
                                  to_screen=True)

    cfg_file_name = sys.argv[1]
    if not os.path.isfile(cfg_file_name):
        msg = "Batch control file " + cfg_file_name + " not found"
        logger.error(msg)
        sys.exit()

    cfg_batch = ConfigObj(cfg_file_name, indent_type="    ", list_values=False,
                          write_empty_values=True)

    main_ui = pfp_batch.Bunch(stop_flag=False, cfg=cfg_batch, mode="batch")

    if cfg_batch["level"] in ["batch_sites"]:
        do_sites_batch(main_ui)
    else:
        msg = " Unrecognised batch type: " + str(cfg_batch["level"])
        logger.error(msg)
