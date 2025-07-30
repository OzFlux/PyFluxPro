# standard modules
import copy
import logging
import os
# 3rd party
import numpy
# PFP modules
from scripts import pfp_utils

pfp_log = os.environ["pfp_log"]
logger = logging.getLogger(pfp_log)

def Standard_deviation_from_variance(ds, Sd_out, Vr_in):
    """
    Purpose:
     Function to convert variance to standard deviation.
    Usage:
     pfp_func_statistics.Standard_deviation_from_variance(ds, Sd_out, Vr_in)
    Author: PRI
    Date: October 2020
    """
    vr_units = {"mg^2/m^6": "mg/m3", "mmol^2/m^6": "mmol/m^3", "g^2/m^6": "g/m^3",
                "degC^2": "degC", "K^2": "K", "m^2/s^2": "m/s"}
    vr = pfp_utils.GetVariable(ds, Vr_in)
    if vr["Attr"]["units"] not in list(vr_units.keys()):
        msg = " Unrecognised units (" + vr["Attr"]["units"] + ") for variable " + Vr_in
        logger.error(msg)
        msg = " Standard deviation not calculated from variance"
        logger.error(msg)
        return 0
    sd = copy.deepcopy(vr)
    sd["Label"] = Sd_out
    sd["Data"] = numpy.ma.sqrt(vr["Data"])
    sd["Attr"]["units"] = vr_units[vr["Attr"]["units"]]
    if "statistic_type" in sd["Attr"]:
        sd["Attr"]["statistic_type"] = "standard_deviation"
    pfp_utils.CreateVariable(ds, sd)
    return 1

def Variance_from_standard_deviation(ds, Vr_out, Sd_in):
    """
    Purpose:
     Function to convert standard deviation to variance.
    Usage:
     pfp_func_statistics.Variance_from_standard_deviation(ds, Vr_out, Sd_in)
    Author: PRI
    Date: October 2020
    """
    sd_units = {"mg/m3": "mg^2/m^6", "mmol/m^3": "mmol^2/m^6", "g/m^3": "g^2/m^6",
                "degC": "degC^2", "K": "K^2", "m/s": "m^2/s^2"}
    sd = pfp_utils.GetVariable(ds, Sd_in)
    if sd["Attr"]["units"] not in list(sd_units.keys()):
        msg = " Unrecognised units (" + sd["Attr"]["units"] + ") for variable " + Sd_in
        logger.error(msg)
        msg = " Variance not calculated from standard deviation"
        logger.error(msg)
        return 0
    vr = copy.deepcopy(sd)
    vr["Label"] = Vr_out
    vr["Data"] = sd["Data"]*sd["Data"]
    vr["Attr"]["units"] = sd_units[sd["Attr"]["units"]]
    if "statistic_type" in vr["Attr"]:
        vr["Attr"]["statistic_type"] = "variance"
    pfp_utils.CreateVariable(ds, vr)
    return 1
