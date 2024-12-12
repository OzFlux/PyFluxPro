# standard modules
import logging
# 3rd party
import dateutil
import numpy
# PFP modules
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def Linear(ds, label_out, label_in, slope, offset, start_date=None, end_date=None):
    # check the input variable exists
    if label_in not in list(ds.root["Variables"].keys()):
        msg = label_in + " not found in data structure, skipping linear ..."
        logger.warning(msg)
        return
    # get the number of records
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    # get the input variable
    var_in = pfp_utils.GetVariable(ds, label_in)
    # create the output variable if it doesn't exist already
    if label_out not in list(ds.root["Variables"].keys()):
        var_out = pfp_utils.CreateEmptyVariable(label_out, nrecs, attr=var_in["Attr"])
    # check the start and end dates
    if start_date is None:
        start_date = var_in["DateTime"][0]
    else:
        try:
            start_date = dateutil.parser.parse(start_date)
        except:
            start_date = var_in["DateTime"][0]
    if end_date is None:
        end_date = var_in["DateTime"][-1]
    else:
        try:
            end_date = dateutil.parser.parse(end_date)
        except:
            end_date = var_in["DateTime"][-1]

    return 1
def Ws_from_Ux_Uy(ds, Ws_out, Ux_in, Uy_in):
    """
    Purpose:
     Function to calculate wind speed from the horizontal components.
    Usage:
     pfp_func_transforms.Ws_from_Ux_Uy(Ux_in, Uy_in)
    Author: PRI
    Date: August 2023
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    attr = {"long_name": "Wind speed", "units": "m/s",
            "statistic_type": "average"}
    Ws = pfp_utils.CreateEmptyVariable(Ws_out, nrecs, attr=attr)
    Ws["Data"] = numpy.ma.sqrt(Ux["Data"]*Ux["Data"] + Uy["Data"]*Uy["Data"])
    Ws["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Ws["Data"]), ones, zeros)
    pfp_utils.CreateVariable(ds, Ws)
    return 1
def Wd_from_Ux_Uy(ds, Wd_out, Ux_in, Uy_in):
    """
    Purpose:
     Function to calculate wind direction from the horizontal components.
    Usage:
     pfp_func_transforms.Wd_from_Ux_Uy(Ux_in, Uy_in)
    Author: PRI
    Date: August 2023
    """
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    attr = {"long_name": "Wind direction", "units": "degrees",
            "statistic_type": "average"}
    Wd = pfp_utils.CreateEmptyVariable(Wd_out, nrecs, attr=attr)
    if ((Ux["Attr"]["instrument"] in ["WindMaster Pro"]) and
        (Uy["Attr"]["instrument"] == Ux["Attr"]["instrument"])):
        Wd_sonic = numpy.ma.mod(360 - numpy.degrees(numpy.ma.arctan2(Uy["Data"], Ux["Data"])), 360)
        Wd["Attr"]["instrument"] = Ux["Attr"]["instrument"]
    elif ((Ux["Attr"]["instrument"] in ["CSAT", "CSAT3A", "CSAT3B"]) and
          (Uy["Attr"]["instrument"] == Ux["Attr"]["instrument"])):
        Wd_sonic = numpy.ma.mod(180 - numpy.degrees(numpy.ma.arctan2(Uy["Data"], Ux["Data"])), 360)
        Wd["Attr"]["instrument"] = Ux["Attr"]["instrument"]
    else:
        msg = " Unrecognised sonic anemometer type (" + Ux["instrument"]
        msg += Uy["instrument"] + ")"
        logger.warning(msg)
    # we want the direction FROM which the wind blows
    Wd["Data"] = numpy.ma.mod(Wd_sonic + 180, 360)
    Wd["Flag"] = numpy.ma.where(numpy.ma.getmaskarray(Wd["Data"]), ones, zeros)
    pfp_utils.CreateVariable(ds, Wd)
    return 1
