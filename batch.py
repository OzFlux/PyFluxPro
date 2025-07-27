# standard modules
import datetime
import logging
import datetime
import logging
import os
import sys
import warnings
# 3rd party modules
from configobj import ConfigObj
# PFP modules
sys.path.append('scripts')
from scripts import pfp_batch
from scripts import pfp_log

warnings.filterwarnings("ignore", category=Warning)

logger = logging.getLogger("pfp_log")

if (__name__ == '__main__'):

    # get the logger
    now = datetime.datetime.now()
    log_file_name = "batch_" + now.strftime("%Y%m%d%H%M") + ".log"
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

    item = pfp_batch.Bunch(stop_flag=False, cfg=cfg_batch, mode="batch")

    if cfg_batch["level"] in ["batch", "batch_levels"]:
        pfp_batch.do_levels_batch(item)
    elif cfg_batch["level"] in ["batch_sites"]:
        pfp_batch.do_sites_batch(item)
    else:
        msg = " Unrecognised batch type: " + str(cfg_batch["level"])
        logger.error(msg)
