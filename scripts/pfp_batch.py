# standard modules
import datetime
import logging
import os
import sys
import traceback
# 3rd party modules
import dateutil
# PFP modules
from scripts import pfp_clim
from scripts import pfp_compliance
from scripts import pfp_cpd_barr
from scripts import pfp_cpd_mchugh
from scripts import pfp_cpd_mcnew
from scripts import pfp_io
from scripts import pfp_levels
from scripts import pfp_mpt
from scripts import pfp_plot
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

class Bunch:
    """
    Constructor class for dummy object with attributes defined by keywords
    when instantiated.
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
def check_file_exists(file_name):
    ok = True
    if not os.path.isfile(file_name):
        msg = file_name + " not found."
        logger.error("")
        logger.error(msg)
        logger.error("")
        ok = False
    return ok
def do_batch_fingerprints(cfg):
    """
    Purpose:
     Plot fingerprints at the end of conatenation, L4 and L5.
    Author: PRI
    Date: Back in the day
    """
    cfg_fp_uri = os.path.join("controlfiles", "standard", "fingerprint.txt")
    cfg_fp = pfp_io.get_controlfilecontents(cfg_fp_uri)
    file_name = pfp_io.get_outfilenamefromcf(cfg)
    file_path = os.path.join(os.path.split(file_name)[0], "")
    plot_path = pfp_utils.get_keyvaluefromcf(cfg, ["Files"], "plot_path", default="plots/")
    cfg_fp["Files"] = {"file_path": file_path, "in_filename": os.path.split(file_name)[1],
                       "plot_path": plot_path}
    cfg_fp["Options"] = {"call_mode": "batch", "show_plots": "No"}
    msg = " Doing fingerprint plots using " + cfg_fp["Files"]["in_filename"]
    logger.info(msg)
    pfp_plot.plot_fingerprint(cfg_fp)
    logger.info("Finished fingerprint plots")
    return
def do_batch_stacked_timeseries(cfg, ds):
    """
    Purpose:
     Plot stacked timeseries for the level being processed.
    Author: PRI
    Date: August 2024
    """
    msg = "Doing stacked timeseries plots"
    logger.info(msg)
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    end = ldt["Data"][-1]
    start = end - dateutil.relativedelta.relativedelta(months=1)
    radn_labels = ["Fsd", "Fsu", "Fld", "Flu", "Fn"]
    flux_labels = ["Fh", "Fe", "Fco2", "Fm"]
    met_labels = ["Precip", "Ta", "RH", "Ws", "Wd", "ps"]
    soil_labels = ["Ts", "Fg", "Sws"]
    plot_labels = radn_labels + flux_labels + met_labels + soil_labels
    cfg["Options"]["plot_stacked_timeseries"] = {"plot_labels": plot_labels,
                                                   "start": start, "end": end,}
    pfp_plot.plot_stacked_timeseries(cfg, ds, start=start, end=end)
    msg = "Finished stacked timeseries plots"
    logger.info(msg)
    return
def do_L1_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L1 processing with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf_l1 = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.l1_update_controlfile(cf_l1, call_mode="batch"):
                msg = "Error occurred updating L1 controlfile"
                logger.error(msg)
                continue
            if not pfp_compliance.check_l1_controlfile(cf_l1):
                msg = "Error occurred checking compliance of L1 controlfile"
                logger.error(msg)
                continue
            ds1 = pfp_levels.l1qc(cf_l1)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l1)
            pfp_io.NetCDFWrite(outfilename, ds1)
            msg = "Finished L1 processing with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during L1 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_L2_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            # get the control file name for this site at this level
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L2 processing with " + cf_file_name[1]
            logger.info(msg)
            # check to see if the control file exists
            if not check_file_exists(cf_level[i]):
                continue
            # get the contents of the control file
            cf_l2 = pfp_io.get_controlfilecontents(cf_level[i])
            # update the L2 control file (legacy)
            if not pfp_compliance.l2_update_controlfile(cf_l2, call_mode="batch"):
                msg = "Error occurred updating L2 controlfile"
                logger.error(msg)
                continue
            # check the contents of the L2 control file
            if not pfp_compliance.check_l2_controlfile(cf_l2):
                msg = "Error occurred checking compliance of L2 controlfile"
                logger.error(msg)
                continue
            # get the input file name
            infilename = pfp_io.get_infilenamefromcf(cf_l2)
            # read the input file
            ds1 = pfp_io.NetCDFRead(infilename)
            # check sonic and IRGA types in L1 netCDF and L2 control file agree
            pfp_compliance.check_l2_options(cf_l2, ds1)
            # do the business
            ds2 = pfp_levels.l2qc(cf_l2, ds1)
            # get the output file name
            outfilename = pfp_io.get_outfilenamefromcf(cf_l2)
            # write the output netCDF file
            pfp_io.NetCDFWrite(outfilename, ds2)
            msg = "Finished L2 processing with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during L2 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_l2_plots_batch(cfg, ds1, ds2):
    if "Plots" in list(cfg.keys()):
        logger.info("Plotting L1 and L2 data")
        for nFig in list(cfg['Plots'].keys()):
            if "(disabled)" in nFig:
                continue
            plt_cf = cfg['Plots'][str(nFig)]
            if 'type' in plt_cf.keys():
                if str(plt_cf['type']).lower() == 'xy':
                    pfp_plot.plotxy(cfg, nFig, plt_cf, ds1, ds2)
                else:
                    pfp_plot.plottimeseries(cfg, nFig, ds1, ds2)
            else:
                pfp_plot.plottimeseries(cfg, nFig, ds1, ds2)
        logger.info("Finished plotting L1 and L2 data")
    return
def do_L3_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L3 processing with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf_l3 = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.l3_update_controlfile(cf_l3, call_mode="batch"):
                msg = "Error occurred updating L3 controlfile"
                logger.error(msg)
                continue
            if not pfp_compliance.check_l3_controlfile(cf_l3):
                msg = "Error occurred checking compliance of L3 controlfile"
                logger.error(msg)
                continue
            infilename = pfp_io.get_infilenamefromcf(cf_l3)
            ds2 = pfp_io.NetCDFRead(infilename)
            ds3 = pfp_levels.l3qc(cf_l3, ds2)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l3)
            pfp_io.NetCDFWrite(outfilename, ds3)
            msg = "Finished L3 processing with " + cf_file_name[1]
            logger.info(msg)
            logger.info("Plotting L3 data")
            do_l3_plots_batch(cf_l3, ds2, ds3)
            logger.info("Finished L3 plotting")
            # plot the L3 fingerprints
            logger.info("Plotting L3 fingerprints")
            do_batch_fingerprints(cf_l3)
            logger.info("Finished L3 fingerprints")
            # do the stacked time series plots
            logger.info("Plotting stacked timeseries")
            do_batch_stacked_timeseries(cf_l3, ds3)
            logger.info("Finished stacked timeseries")
            logger.info("")
        except Exception:
            msg = "Error occurred during L3 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_l3_plots_batch(cfg, ds2, ds3):
    if "Plots" in list(cfg.keys()):
        for nFig in list(cfg['Plots'].keys()):
            if "(disabled)" in nFig:
                continue
            plt_cf = cfg['Plots'][str(nFig)]
            if 'type' in plt_cf.keys():
                if str(plt_cf['type']).lower() == 'xy':
                    pfp_plot.plotxy(cfg, nFig, plt_cf, ds2, ds3)
                else:
                    pfp_plot.plottimeseries(cfg, nFig, ds2, ds3)
            else:
                pfp_plot.plottimeseries(cfg, nFig, ds2, ds3)
    return
def do_reddyproc_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting REddyProc output with " + cf_file_name[1]
        logger.info(msg)
        if not check_file_exists(cf_level[i]):
            continue
        cf = pfp_io.get_controlfilecontents(cf_level[i])
        pfp_io.write_tsv_reddyproc(cf)
        msg = "Finished REddyProc output with " + cf_file_name[1]
        logger.info(msg)
        logger.info("")
    return ok
def do_oneflux_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting ONEFlux output with " + cf_file_name[1]
        logger.info(msg)
        if not check_file_exists(cf_level[i]):
            continue
        cf = pfp_io.get_controlfilecontents(cf_level[i])
        pfp_io.write_csv_oneflux(cf)
        msg = "Finished ONEFlux output with " + cf_file_name[1]
        logger.info(msg)
        logger.info("")
    return ok
def do_concatenate_batch(main_ui, cf_level):
    ok = True
    sites = sorted(list(cf_level.keys()), key=int)
    for i in sites:
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            if not check_file_exists(cf_level[i]):
                continue
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting concatenation with " + cf_file_name[1]
            logger.info(msg)
            cf_cc = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.concatenate_update_controlfile(cf_cc):
                msg = "Error occurred updating concatenation controlfile"
                logger.error(msg)
                continue
            info = pfp_compliance.ParseConcatenateControlFile(cf_cc)
            if not info["NetCDFConcatenate"]["OK"]:
                msg = " Error occurred parsing the control file " + cf_file_name[1]
                logger.error(msg)
                continue
            pfp_io.NetCDFConcatenate(info)
            msg = "Finished concatenation with " + cf_file_name[1]
            logger.info(msg)
            # and then plot the fingerprints for the concatenated files
            do_batch_fingerprints(cf_cc)
            logger.info("")
        except Exception:
            msg = "Error occurred during concatenation with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_climatology_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        if not check_file_exists(cf_level[i]):
            continue
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting climatology with " + cf_file_name[1]
            logger.info(msg)
            cf_ct = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.climatology_update_controlfile(cf_ct):
                msg = "Error occurred updating climatology controlfile"
                logger.error(msg)
                continue
            pfp_clim.climatology(cf_ct)
            msg = "Finished climatology with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during climatology with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_cpd_barr_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting CPD (Barr) with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.cpd_barr_update_controlfile(cf, call_mode="batch"):
                msg = "Error occurred updating CPD (Barr) controlfile"
                logger.error(msg)
                continue
            pfp_cpd_barr.cpd_barr_main(cf)
            msg = "Finished CPD (Barr) with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during CPD (Barr) with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_cpd_mchugh_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting CPD (McHugh) with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.cpd_mchugh_update_controlfile(cf, call_mode="batch"):
                msg = "Error occurred updating CPD (McHugh) controlfile"
                logger.error(msg)
                continue
            pfp_cpd_mchugh.cpd_mchugh_main(cf)
            msg = "Finished CPD (McHugh) with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during CPD (McHugh) with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_cpd_mcnew_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting CPD (McNew) with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.cpd_mcnew_update_controlfile(cf, call_mode="batch"):
                msg = "Error occurred updating CPD (McNew) controlfile"
                logger.error(msg)
                continue
            pfp_cpd_mcnew.cpd_mcnew_main(cf)
            msg = "Finished CPD (McNew) with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during CPD (McNew) with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_mpt_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting MPT with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.mpt_update_controlfile(cf, call_mode="batch"):
                msg = "Error occurred updating MPT controlfile"
                logger.error(msg)
                continue
            pfp_mpt.mpt_main(cf)
            msg = "Finished MPT with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during MPT with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_L4_batch(main_ui, cf_level):
    ok = True
    sites = sorted(list(cf_level.keys()), key=int)
    for i in sites:
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L4 processing with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf_l4 = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.l4_update_controlfile(cf_l4, call_mode="batch"):
                msg = "Error occurred updating L4 controlfile"
                logger.error(msg)
                continue
            infilename = pfp_io.get_infilenamefromcf(cf_l4)
            ds3 = pfp_io.NetCDFRead(infilename)
            ds4 = pfp_levels.l4qc(None, cf_l4, ds3)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l4)
            pfp_io.NetCDFWrite(outfilename, ds4)
            msg = "Finished L4 processing with " + cf_file_name[1]
            logger.info(msg)
            # plot the L4 fingerprints
            do_batch_fingerprints(cf_l4)
            logger.info("")
        except Exception:
            msg = "Error occurred during L4 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_L5_batch(main_ui, cf_level):
    ok = True
    sites = sorted(list(cf_level.keys()), key=int)
    for i in sites:
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L5 processing with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf_l5 = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.l5_update_controlfile(cf_l5, call_mode="batch"):
                msg = "Error occurred updating L5 controlfile"
                logger.error(msg)
                continue
            if not pfp_compliance.check_l5_controlfile(cf_l5):
                msg = "Error occurred checking compliance of L5 controlfile"
                logger.error(msg)
                continue
            infilename = pfp_io.get_infilenamefromcf(cf_l5)
            ds4 = pfp_io.NetCDFRead(infilename)
            ds5 = pfp_levels.l5qc(None, cf_l5, ds4)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l5)
            pfp_io.NetCDFWrite(outfilename, ds5)
            msg = "Finished L5 processing with " + cf_file_name[1]
            logger.info(msg)
            # plot the L5 fingerprints
            do_batch_fingerprints(cf_l5)
            logger.info("")
        except Exception:
            msg = "Error occurred during L5 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_L6_batch(main_ui, cf_level):
    ok = True
    for i in list(cf_level.keys()):
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        try:
            cf_file_name = os.path.split(cf_level[i])
            msg = "Starting L6 processing with " + cf_file_name[1]
            logger.info(msg)
            if not check_file_exists(cf_level[i]):
                continue
            cf_l6 = pfp_io.get_controlfilecontents(cf_level[i])
            if not pfp_compliance.l6_update_controlfile(cf_l6, call_mode="batch"):
                msg = "Error occurred updating L6 controlfile"
                logger.error(msg)
                continue
            if not pfp_compliance.check_l6_controlfile(cf_l6):
                msg = "Error occurred checking compliance of L6 controlfile"
                logger.error(msg)
                continue
            infilename = pfp_io.get_infilenamefromcf(cf_l6)
            ds5 = pfp_io.NetCDFRead(infilename)
            ds6 = pfp_levels.l6qc(None, cf_l6, ds5)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l6)
            pfp_io.NetCDFWrite(outfilename, ds6)
            msg = "Finished L6 processing with " + cf_file_name[1]
            logger.info(msg)
            # do the CF compliance check
            #do_batch_cfcheck(cf_l6)
            logger.info("")
        except Exception:
            msg = "Error occurred during L6 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return ok
def do_levels_batch(main_ui):
    logger = logging.getLogger("pfp_log")
    if main_ui.mode == "interactive":
        tab_index_running = main_ui.tabs.tab_index_running
        cf_batch = main_ui.tabs.tab_dict[tab_index_running].get_data_from_model()
    elif main_ui.mode == "batch":
        cf_batch = main_ui.cfg
    else:
        msg = "Unrecognised option for mode (" + main_ui.mode + ")"
        logger.error(msg)
        raise RuntimeError
    start = datetime.datetime.now()
    msg = "Started batch processing at " + start.strftime("%Y%m%d%H%M")
    logger.info(msg)
    if "Options" in cf_batch:
        if "levels" in cf_batch["Options"]:
            levels = pfp_utils.string_to_list(cf_batch["Options"]["levels"])
        else:
            msg = "No 'levels' entry found in [Options] section"
            logger.error(msg)
            sys.exit()
    else:
        msg = "No [Options] section in control file"
        logger.error(msg)
        sys.exit()
    processing_levels = ["l1", "l2", "l3",
                         "nc2csv_oneflux",
                         "concatenate", "climatology",
                         "cpd_barr", "cpd_mchugh", "cpd_mcnew", "mpt",
                         "l4", "l5", "l6"]
    for level in levels:
        # check the stop flag
        if main_ui.stop_flag:
            # break out of the loop if user requested stop
            break
        if level.lower() not in processing_levels:
            msg = "Unrecognised level " + level
            logger.warning(msg)
            continue
        if level.lower() == "l1":
            # L1 processing
            if not do_L1_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "l2":
            # L2 processing
            if not do_L2_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "l3":
            # L3 processing
            if not do_L3_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "nc2csv_oneflux":
            # convert netCDF files to ONEFlux CSV files
            if not do_oneflux_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "reddyproc":
            # convert netCDF files to REddyProc CSV files
            if not do_reddyproc_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "concatenate":
            # concatenate netCDF files
            if not do_concatenate_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "climatology":
            # climatology
            if not do_climatology_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "cpd_barr":
            # ustar threshold from change point detection
            if not do_cpd_barr_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "cpd_mchugh":
            # ustar threshold from change point detection
            if not do_cpd_mchugh_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "cpd_mcnew":
            # ustar threshold from change point detection
            if not do_cpd_mcnew_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "mpt":
            # ustar threshold from change point detection
            if not do_mpt_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "l4":
            # L4 processing
            if not do_L4_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "l5":
            # L5 processing
            if not do_L5_batch(main_ui, cf_batch["Levels"][level]):
                break
        elif level.lower() == "l6":
            # L6 processing
            if not do_L6_batch(main_ui, cf_batch["Levels"][level]):
                break
    end = datetime.datetime.now()
    msg = " Finished batch processing at " + end.strftime("%Y%m%d%H%M")
    logger.info(msg)
    return
def do_sites_batch(item):
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
    # loop over the sites
    for site in list(item.cfg["Sites"].keys()):
        # get the control files for this site
        cfg_site = item.cfg["Sites"][site]
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
                do_L1_batch(item, cf_level)
            elif level.lower() == "l2":
                # L2 processing
                do_L2_batch(item, cf_level)
            elif level.lower() == "l3":
                # L3 processing
                do_L3_batch(item, cf_level)
            elif level.lower() == "concatenate":
                # concatenate netCDF files
                do_concatenate_batch(item, cf_level)
            elif level.lower() == "climatology":
                # climatology
                do_climatology_batch(item, cf_level)
            elif level.lower() == "cpd_barr":
                # ustar threshold from change point detection
                do_cpd_barr_batch(item, cf_level)
            elif level.lower() == "cpd_mchugh":
                # ustar threshold from change point detection
                do_cpd_mchugh_batch(item, cf_level)
            elif level.lower() == "cpd_mcnew":
                # ustar threshold from change point detection
                do_cpd_mcnew_batch(item, cf_level)
            elif level.lower() == "mpt":
                # ustar threshold from change point detection
                do_mpt_batch(item, cf_level)
            elif level.lower() == "l4":
                # L4 processing
                do_L4_batch(item, cf_level)
            elif level.lower() == "l5":
                # L5 processing
                do_L5_batch(item, cf_level)
            elif level.lower() == "l6":
                # L6 processing
                do_L6_batch(item, cf_level)
            else:
                msg = " Unrecognised batch processing level " + str(level)
                logger.error(msg)
