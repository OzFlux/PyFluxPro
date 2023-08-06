# standard modules
import logging
# 3rd party
import numpy
# PFP modules
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def Ws_from_Ux_Uy(ds, Ws_out, Ux_in, Uy_in):
    """
    Purpose:
     Function to calculate wind speed from the horizontal components.
    Usage:
     pfp_func_transforms.Ws_from_Ux_Uy(Ux_in, Uy_in)
    Author: PRI
    Date: August 2023
    """
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    Ws, _ = pfp_utils.convert_UVtoWSWD(Ux, Uy)
    Ws["Label"] = Ws_out
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
    Ux = pfp_utils.GetVariable(ds, Ux_in)
    Uy = pfp_utils.GetVariable(ds, Uy_in)
    _, Wd = pfp_utils.convert_UVtoWSWD(Ux, Uy)
    Wd["Label"] = Wd_out
    pfp_utils.CreateVariable(ds, Wd)
    return 1
