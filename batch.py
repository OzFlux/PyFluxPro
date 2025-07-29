# standard modules
import datetime
import os
import sys
import warnings
# 3rd party modules
from configobj import ConfigObj
# PFP modules
#sys.path.append('scripts')
#from scripts import pfp_batch
#from scripts import pfp_log

warnings.filterwarnings("ignore", category=Warning)

if (__name__ == '__main__'):

    cfg_file_url = sys.argv[1]
    if not os.path.isfile(cfg_file_url):
        msg = "Batch control file " + cfg_file_url + " not found"
        print(msg)
        sys.exit()

    # get the logger
    now = datetime.datetime.now()
    cfg_base_name, _= os.path.splitext(os.path.basename(cfg_file_url))
    log_file_base = cfg_base_name + "_" + now.strftime("%Y%m%dT%H%M%S%f")
    os.environ["pfp_log"] = log_file_base
    log_file_name = log_file_base + ".log"
    log_file_url = os.path.join("logfiles", log_file_name)

    from scripts import pfp_batch
    from scripts import pfp_log

    logger = pfp_log.CreateLogger(log_file_base, log_file_name=log_file_url, to_screen=True)

    cfg_batch = ConfigObj(cfg_file_url, indent_type="    ", list_values=False, write_empty_values=True)

    item = pfp_batch.Bunch(stop_flag=False, cfg=cfg_batch, mode="batch")

    if cfg_batch["level"] in ["batch", "batch_levels"]:
        pfp_batch.do_levels_batch(item)
    elif cfg_batch["level"] in ["batch_sites"]:
        pfp_batch.do_sites_batch(item)
    else:
        msg = " Unrecognised batch type: " + str(cfg_batch["level"])
        logger.error(msg)
