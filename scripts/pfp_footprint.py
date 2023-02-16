import datetime
import logging
import math
import os
import copy
import warnings
# 3rd party
import matplotlib
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    matplotlib.rcParams['toolbar'] = 'toolmanager'
import matplotlib.dates as mdt
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backend_tools import ToolToggleBase
import numpy
from scipy import stats
import statsmodels.api as sm
# PFP modules
from scripts import constants as c
from scripts import pfp_classes
from scripts import pfp_ck
from scripts import pfp_io
from scripts import pfp_utils
from scripts import pfp_ts
from scripts import meteorologicalfunctions as pfp_mf

"""
Derive a flux footprint estimate based on the paper from Kormann and Meixner (2001). The equations solved
follow the work from Neftel et al. (2008) and the ART footprint tool. Currently no sections are 
defined to estimate the proportion of footprint source from this defined areas as in Neftel et al. (2008).
See Kormann, R. and Meixner, F.X., 2001: An analytical footprint model for non-neutral stratification.
Boundary-Layer Meteorology 99: 207. https://doi.org/10.1023/A:1018991015119 and 
Neftel, A., Spirig, C., Ammann, C., 2008: Application and test of a simple tool for operational footprint 
evaluations. Environmental Pollution 152, 644-652. for details.
The layout of this python script follows the footprint climatology layout by Kljun et al., 2015.
contact: cacilia.ewenz@internode.on.net

This function calculates footprints within a fixed physical domain for a series of
time steps, rotates footprints into the corresponding wind direction and aggregates
all footprints to a footprint climatology. The percentage of source area is
calculated for the footprint climatology.


FKM Input
    All vectors need to be of equal length (one value for each time step)
    zm       = Measurement height above displacement height (i.e. z-d) [m]
               usually a scalar, but can also be a vector 
    z0       = Roughness length [m]
               usually a scalar, but can also be a vector 
    umean    = Vector of mean wind speed at zm [ms-1] - enter [None] if not known 
               Either z0 or umean is required. If both are given,
               z0 is selected to calculate the footprint
    ol       = Vector of Obukhov length [m]
    sigmav   = Vector of standard deviation of lateral velocity fluctuations [ms-1]
               if not available it is calculated as 0.5*Ws (wind speed)
    ustar    = Vector of friction velocity [ms-1]
    wind_dir = Vector of wind direction in degrees (of 360) for rotation of the footprint

    Optional input:
    domain       = Domain size as an array of [xmin xmax ymin ymax] [m].
                   Footprint will be calculated for a measurement at [0 0 zm] m
                   Default is smallest area including the r% footprint or [-1000 1000 -1000 1000]m,
                   whichever smallest (80% footprint if r not given).
    dx, dy       = Cell size of domain [m]
                   Small dx, dy results in higher spatial resolution and higher computing time
                   Default is dx = dy = 2 m. If only dx is given, dx=dy.
    nx, ny       = Two integer scalars defining the number of grid elements in x and y
                   Large nx/ny result in higher spatial resolution and higher computing time
                   Default is nx = ny = 1000. If only nx is given, nx=ny.
                   If both dx/dy and nx/ny are given, dx/dy is given priority if the domain is also specified.
    pulse        = Display progress of footprint calculations every pulse-the footprint (e.g., "100")
    verbosity    = Level of verbosity at run time: 0 = completely silent, 1 = notify only of fatal errors,
                   2 = all notifications
FKM output
    FKM      = Structure array with footprint climatology data for measurement at [0 0 zm] m
    x_2d     = x-grid of 2-dimensional footprint [m]
    y_2d     = y-grid of 2-dimensional footprint [m]
    fclim_2d = Normalised footprint function values of footprint climatology [m-2]
    n        = Number of footprints calculated and included in footprint climatology
    flag_err = 0 if no error, 1 in case of error, 2 if not all contour plots (rs%) within specified domain

Footprint calculation following the Kormann-Meixner Approach
1st Version by Marx Stampfli, stampfli mathematics, Bern, Switzerland
Revision December 2006, C. Spirig and A. Neftel, Agroscope, Zurich, Switzerland
            ported into python CM Ewenz, Adelaide Nov 2015 to May 2016
Array of lower and upper ranges for following parameters
             (empty,Xcoord,Ycoord,U_star,     LM,Std_v,Wdir,   zm,U_meas,empty)
lrange = Array(0,    -5000, -5000,  0.01,-999999,    0,   0,    0,     0,    0)
urange = Array(0,     5000,  5000,     5, 999999,   20, 360, 1000,    30,    0)

Created: May 2018 Cacilia Ewenz
version: 0.1
last change: 08/06/2018 Cacilia Ewenz
Copyright (C) 2018, Cacilia Ewenz
"""
logger = logging.getLogger("pfp_log")

def calculate_footprint(cf):
    """
    Purpose:
    Usage:
    Side effects:
    Author: PRI/CE
    Date: Feb 2023
    """
    # read the controlfile
    fpinfo = calculate_footprint_parse_controlfile(cf)
    # read the netCDF file to a data structure
    ds = pfp_io.NetCDFRead(fpinfo["in_filename"])
    # calculate Kormann-Meixner footprint
    ds_fkm = pfp_classes.Datastructure()
    # copy global attributes to new dictionary
    ds_fkm = copy.deepcopy(ds)
    # check the required variables are in the data structure
    if not calculate_footprint_check_variables(ds_fkm, fpinfo):
        return
    # update the info dictionary
    fpinfo["start"] = ds_fkm.root["Attributes"]["time_coverage_start"]
    fpinfo["end"] = ds_fkm.root["Attributes"]["time_coverage_end"]
    fpinfo["site_name"] = ds_fkm.root["Attributes"]["site_name"]
    # (1) Evaluate input data and apply constraints
    check_fkm_inputs(ds_fkm,fpinfo)
    # (2) get a mask for all input variables
    composite_mask(ds_fkm, fpinfo)
    # (3) calculate Kormann and Meixner, 2001 footprint
    kormann_meixner(ds_fkm, fpinfo)
    # plot all the data
    #plot_footprint_all(ds, fpinfo)
    # plot the seasonal footprint
    #plot_footprint_seasonal(ds, fpinfo)
    return

def calculate_footprint_check_variables(ds, fpinfo):
    """ Check the variables named in the control file exist in the data structure."""
    ok = True
    labels = list(ds.root["Variables"].keys())
    for item in ["Fsd", "Wd", "Ws", "ustar", "Vsd", "ol"]:
        name = fpinfo[item]
        if name not in labels:
            msg = name + "(" + item + ") not found in " + fpinfo["in_filename"]
            logger.error(msg)
            ok = False
            if item == "ol":
                msg = name + "(" + item + ") will be calculated."
                logger.error(msg)
                pfp_ts.CalculateMoninObukhovLength(ds)
                ok = True
            if item == "Vsd":
                msg = name + "(" + item + ") will be calculated."
                logger.error(msg)
                crosswind_std(ds, fpinfo)
                ok = True
    return ok

def calculate_footprint_parse_controlfile(cf):
    """
    Purpose:
     Parse the footprint control file and return the required information in a
     dictionary.
    Author: PRI/CE
    Date: Feb 2023
    """
    fpinfo = {}
    # get the day/night separator Fsd name and its threshold
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "Fsd_threshold", default="10")
    fpinfo["Fsd_threshold"] = int(opt)
    fpinfo["Fsd"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Fsd"], "name", default="Fsd")
    # get the wind direction and wind speed variable names
    fpinfo["Wd"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Wd"], "name", default="Wd")
    fpinfo["Ws"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Ws"], "name", default="Ws")
    # friction velocity
    fpinfo["ustar"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "ustar"], "name", default="ustar")
    # Monin-Obukhov-length if available
    fpinfo["ol"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "ol"], "name", default="L")
    # cross wind speed standard deviation
    fpinfo["Vsd"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Vsd"], "name", default="Vsd")
    # plot path
    fpinfo["plot_path"] = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "plot_path", default="plots/")
    if not os.path.exists(fpinfo["plot_path"]):
        os.makedirs(fpinfo["plot_path"])
    return fpinfo

def check_fkm_inputs(ds,fpinfo):
    """
    Purpose:
     Check passed values for physical plausibility and consistency and reduce the timeseries of only valid input data.
    Author: CE
    Date: Feb 2023
    """
    logger.info('Check range for input variables')
    
    # get the required variables
    ol = pfp_utils.GetVariable(ds, fpinfo["ol"])
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])
    
    # get the instrument height above the displacement height (2/3*canopy_height)
    th = float(pfp_utils.strip_non_numeric(ds.root["Attributes"]["tower_height"]))
    dh = float(pfp_utils.strip_non_numeric(ds.root["Attributes"]["canopy_height"]))
    fpinfo["zm"] = th - (2/3)*dh
    # get the roughness length and z/L
    z0 = z0calc(ds, fpinfo)
    
    # create a variable dictionary for zmol
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    zmol = pfp_utils.CreateEmptyVariable("zmol", nrecs, datetime=ldt["Data"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    zmol["Data"] = fpinfo["zm"]/ol["Data"]
    # get the QC flag
    zmol["Flag"] = numpy.where(numpy.ma.getmaskarray(zmol["Data"]) == True, ones, zeros)
    # update the variable attributes
    zmol["Data"]["Attr"]["units"] = "1"
    zmol["Data"]["Attr"]["long_name"] = "zm over L"
    # put the cross wind standard deviation variable in the data structure
    pfp_utils.CreateVariable(ds, zmol)
    
    # start checking constraints for footprint calculation
    if fpinfo["zm"] <= 0.:
        raise_fkm_exception(2, verbosity, counter, msgstring)
        return False
    #if z0 <= 0.:
    #    raise_fkm_exception(3, verbosity, counter, msgstring)
    #    return False
    z0["Data"] = numpy.ma.masked_where(ds.root["Variables"]["z0"] < 0.0, z0["Data"])
    #if float(zm)/ol < -3:
    #    raise_fkm_exception(7, verbosity, counter, msgstring)
    #    return False
    # create zm/L
    zmol["Data"] = numpy.ma.masked_where(zmol["Data"] < -3.0, zmol["Data"])
    #if float(zm)/ol > 3:
    #    raise_fkm_exception(7, verbosity, counter, msgstring)
    #    return False
    zmol["Data"] = numpy.ma.masked_where(zmol["Data"] > 3.0, zmol["Data"])
    #if sigmav <= 0:
    #    raise_fkm_exception(8, verbosity, counter, msgstring)
    #    return False
    Vsd["Data"] = numpy.ma.masked_where(Vsd["Data"] <= 0.0, Vsd["Data"])
    #if ustar <= 0.1:
    #    raise_fkm_exception(9, verbosity, counter, msgstring)
    #    return False
    ustar["Data"] = numpy.ma.masked_where(ustar["Data"] <= 0.1, ustar["Data"])
    #if umean <= 0.0:
    #    raise_fkm_exception(21, verbosity, counter, msgstring)
    #    return False
    Ws["Data"] = numpy.ma.masked_where(Ws["Data"] <= 0.1, Ws["Data"])
    #if wind_dir > 360:
    #    raise_fkm_exception(10, verbosity, counter, msgstring)
    #    return False
    Wd["Data"] = numpy.ma.masked_where(Wd["Data"] > 360.0, Wd["Data"])
    #if wind_dir < 0:
    #    raise_fkm_exception(10, verbosity, counter, msgstring)
    #    return False
    Wd["Data"] = numpy.ma.masked_where(Wd["Data"] < 0.0, Wd["Data"])
    return True

def composite_mask(ds, fpinfo):
    """
    Purpose:
     Calculate and apply a composite mask for all Kormann-Maixner input variables.
    Author: CE
    Date: Feb 2023
    """
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    # get the required variables
    z0 = pfp_utils.GetVariable(ds, fpinfo["z0"])
    zmol = pfp_utils.GetVariable(ds, fpinfo["zmol"])
    ol = pfp_utils.GetVariable(ds, fpinfo["ol"])
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])

    # get a composite mask over all variables needed for Kormann-Meixner 2001 footprint
    mask_all = numpy.ma.getmaskarray(ds.root["Variables"]["ustar"])
    for item in ["Ws", "ol", "Vsd", "Wd", "z0", "zmol"]:
        mask_item = numpy.ma.getmaskarray(ds.root["Variables"][item])
        mask_all = numpy.ma.mask_or(mask_all, mask_item)
    
    # and then apply the composite mask to all variables and remove masked elements
    Ws["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, Ws["Data"])))
    ol["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ol["Data"])))
    zmol["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, zmol["Data"])))
    Vsd["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, Vsd["Data"])))
    ustar["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ustar["Data"])))
    Wd["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, Wd["Data"])))
    z0["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, z0["Data"])))
    ldt["Data"] = list(numpy.ma.compressed(numpy.ma.masked_where(mask_all == True, ldt["Data"])))
    if len(ustar["Data"]) == 0:
        msg = "No footprint input data for "+str(ldt[0])+" to "+str(ldt[-1])
        logger.warning(msg)
    return

def kormann_meixner(ds, fpinfo):
    """
    ustar, sigmav, ol, wind_dir, zm, z0, umean
    Purpose:
     Calculate Kormann-Maixner footprint according to Neftel et al., 2008
     ART Footprint Model.
    Author: CE
    Date: Feb 2023
    """
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    nrecs = ldt["Data"]
    # get the required variables
    z0 = pfp_utils.GetVariable(ds, fpinfo["z0"])
    ol = pfp_utils.GetVariable(ds, fpinfo["ol"])
    zmol = pfp_utils.GetVariable(ds, fpinfo["zmol"])
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])

    zm = fpinfo["zm"]
    # a) create positive, zero and negative LM masks
    LMp = numpy.ma.masked_where(ol["Data"] <  float(0),ol["Data"])
    LMz = numpy.ma.masked_where(ol["Data"] == float(0),ol["Data"])
    LMn = numpy.ma.masked_where(ol["Data"] > float(0),ol["Data"])
    # phi, u, m, n
    p_phi = 1.0+5.0*zmol["Data"]
    p_u   = (ustar["Data"] / c.k) * (numpy.log(zm / z0) + 5.0 * zmol["Data"])
    p_m   = (1.0 + 5.0 * (zmol["Data"])) / (numpy.log(zm / z0["Data"]) + 5.0 * zmol["Data"])
    p_n   = 1.0 / (1.0 + 5.0 * zmol["Data"])
    z_phi = 1.0
    z_u   = ustar["Data"] / c.k * (numpy.log(zm / z0["Data"]))
    z_m   = ustar["Data"] / c.k / z_u
    z_n   = 1.0
    zeta = (1.0-16.0*zmol["Data"])**0.25
    psi = -2.0*numpy.log(0.5*(1.0+zeta))-numpy.log(0.5*(1.0+zeta*zeta))+2.0*numpy.arctan(zeta)-0.5*c.Pi
    n_phi = 1.0/(zeta**2)
    n_u   = ustar["Data"]/c.k*(numpy.log(zm/z0["Data"])+psi)
    n_m   = ustar["Data"]/c.k/zeta/n_u
    n_n   = (1.0-24.0*zmol["Data"])/(1.0-16.0*zmol["Data"])
    # mask the arrays
    p_phi = numpy.ma.masked_where(LMp==True,p_phi)
    p_u   = numpy.ma.masked_where(LMp==True,p_u)
    p_m   = numpy.ma.masked_where(LMp==True,p_m)
    p_n   = numpy.ma.masked_where(LMp==True,p_n)
    z_phi = numpy.ma.masked_where(LMz==True,z_phi)
    z_u   = numpy.ma.masked_where(LMz==True,z_u)
    z_m   = numpy.ma.masked_where(LMz==True,z_m)
    z_n   = numpy.ma.masked_where(LMz==True,z_n)
    n_phi = numpy.ma.masked_where(LMn==True,n_phi)
    n_u   = numpy.ma.masked_where(LMn==True,n_u)
    n_m   = numpy.ma.masked_where(LMn==True,n_m)
    n_n   = numpy.ma.masked_where(LMn==True,n_n)
    # fill positive, zero and negative masks
    p_phi = numpy.ma.filled(p_phi,float(0))
    p_u   = numpy.ma.filled(p_u,float(0))
    p_m   = numpy.ma.filled(p_m,float(0))
    p_n   = numpy.ma.filled(p_n,float(0))
    z_phi = numpy.ma.filled(z_phi,float(0))
    z_u   = numpy.ma.filled(z_u,float(0))
    z_m   = numpy.ma.filled(z_m,float(0))
    z_n   = numpy.ma.filled(z_n,float(0))
    n_phi = numpy.ma.filled(n_phi,float(0))
    n_u   = numpy.ma.filled(n_u,float(0))
    n_m   = numpy.ma.filled(n_m,float(0))
    n_n   = numpy.ma.filled(n_n,float(0))
    # add the arrays together to get a full array
    phi = p_phi + z_phi + n_phi
    u   = p_u + z_u + n_u
    m   = p_m + z_m + n_m
    n   = p_n + z_n + n_n
    # --- calculate Korman-Meixner footprint
    # r, mu
    r = 2 + m - n
    mu = (1 + m) / (2 + m - n)
    # U=Umaj, Kmaj, xi
    Umaj = u / (zm**m)
    Kmaj = c.k * ustar * zm / phi / zm**n
    # Kmaj corresponds to kappa in KM 2001
    xi = Umaj * (zm **r) / (r * r * Kmaj)
    # xPhiMax is the (x-)position of the maximum of phi
    xPhiMax = r * xi / (2 * r + 1)
    # determine GammaProxmu and GammaProx1r depending on mu and 1/r
    GammaProxmu = numpy.zeros_like(mu)
    GammaProxmu = numpy.ma.masked_where(mu == 0.0, GammaProxmu)
    GammaProxmu = (1.0 / mu) + 0.1002588 * numpy.exp(mu) - 0.493536 + 0.3066 * mu - 0.09 * (mu**2)
    GammaProxmu   = numpy.ma.filled(GammaProxmu,float(0))
    GammaProx1r = numpy.zeros_like(r)
    GammaProx1r = numpy.ma.masked_where(1/r == 0.0, GammaProx1r)
    GammaProx1r = r + 0.1002588 * numpy.exp(1/r) - 0.493536 + 0.3066 * (1/r) - 0.09 * ((1/r)**2)
    GammaProx1r   = numpy.ma.filled(GammaProx1r,float(0))
    # Kormann-Meixner parameters A-E
    A = pfp_utils.CreateEmptyVariable("A", nrecs, datetime=ldt["Data"])
    B = pfp_utils.CreateEmptyVariable("B", nrecs, datetime=ldt["Data"])
    C = pfp_utils.CreateEmptyVariable("C", nrecs, datetime=ldt["Data"])
    D = pfp_utils.CreateEmptyVariable("D", nrecs, datetime=ldt["Data"])
    E = pfp_utils.CreateEmptyVariable("E", nrecs, datetime=ldt["Data"])
    KM_p01a = pfp_utils.CreateEmptyVariable("KM_p01a", nrecs, datetime=ldt["Data"])
    KM_p01b = pfp_utils.CreateEmptyVariable("KM_p01b", nrecs, datetime=ldt["Data"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    A["Data"] = 1 + mu
    B["Data"] = Umaj * (zm**r) / r / r / Kmaj
    C["Data"] = (B**mu) /GammaProxmu
    D["Data"] = Vsd["Data"] * GammaProx1r / GammaProxmu / ((r * r * Kmaj / Umaj)**(m / r)) / Umaj
    E["Data"] = (r - m) / r
    # Ellipse parameters: Calculate the near and far distance for 1% of the max level
    philevel = 0.01
    x0 = xPhiMax / 3
    p = philevel * xPhiMax
    for _ in [1,2,3,4,5]:
        Fm = C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"]) - p
        dFm = (C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"])) / x0**2 * -(A["Data"] * x0 - B["Data"]) - p
        x0 = x0 - Fm / dFm
    KM_p01a["Data"] = x0
    x0 = xPhiMax * 3
    p = philevel * xPhiMax
    for _ in [1,2,3,4,5]:
        Fm = C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"]) - p
        dFm = (C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"])) / x0**2 * -(A["Data"] * x0 - B["Data"]) - p
        x0 = x0 - Fm / dFm
    KM_p01b["Data"] = x0
    # get the QC flag
    A["Flag"] = numpy.where(numpy.ma.getmaskarray(A["Data"]) == True, ones, zeros)
    B["Flag"] = numpy.where(numpy.ma.getmaskarray(B["Data"]) == True, ones, zeros)
    C["Flag"] = numpy.where(numpy.ma.getmaskarray(C["Data"]) == True, ones, zeros)
    D["Flag"] = numpy.where(numpy.ma.getmaskarray(D["Data"]) == True, ones, zeros)
    E["Flag"] = numpy.where(numpy.ma.getmaskarray(E["Data"]) == True, ones, zeros)
    KM_p01a["Flag"] = numpy.where(numpy.ma.getmaskarray(KM_p01a["Data"]) == True, ones, zeros)
    KM_p01b["Flag"] = numpy.where(numpy.ma.getmaskarray(KM_p01b["Data"]) == True, ones, zeros)
    # update the variable attributes
    A["Data"]["Attr"]["units"] = "1"
    B["Data"]["Attr"]["units"] = "1"
    C["Data"]["Attr"]["units"] = "1"
    D["Data"]["Attr"]["units"] = "1"
    E["Data"]["Attr"]["units"] = "1"
    KM_p01a["Data"]["Attr"]["units"] = "m"
    KM_p01b["Data"]["Attr"]["units"] = "m"
    A["Data"]["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    B["Data"]["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    C["Data"]["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    D["Data"]["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    E["Data"]["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    KM_p01a["Data"]["Attr"]["long_name"] = "Kormann-Meixner ellipse 1% value close distance from the tower"
    KM_p01b["Data"]["Attr"]["long_name"] = "Kormann-Meixner ellipse 1% value far distance from the tower"
    pfp_utils.CreateVariable(ds, A)
    pfp_utils.CreateVariable(ds, B)
    pfp_utils.CreateVariable(ds, C)
    pfp_utils.CreateVariable(ds, D)
    pfp_utils.CreateVariable(ds, E)
    pfp_utils.CreateVariable(ds, KM_p01a)
    pfp_utils.CreateVariable(ds, KM_p01b)
    return

def raise_fkm_exception(code, verbosity, counter, msgstring):
    '''Raise exception or prints message according to specified code'''

    icode = int(code)
    if icode > 0:
        if counter[icode] == None:
            counter[icode] = 1
        else:
            counter[icode] = counter[icode] + 1
        ex = [it for it in exceptions if it['code'] == code][0]
        msgstring[icode] = ex['type'] + '(' + str(ex['code']).zfill(4) + '):\n '+ ex['msg']
        #if verbosity > 0: print('')
        if ex['type'] == exTypes['fatal']:
            if verbosity > 0:
                msgstring[icode] = msgstring[icode] + '\n FKM_fixed_domain execution aborted.'
            else:
                msgstring[icode] = ''
            raise Exception(msgstring[icode])
        elif ex['type'] == exTypes['alert']:
            msgstring[icode] = msgstring[icode] #+ '\n Execution continues.'
            if verbosity > 1: 
                pass
        elif ex['type'] == exTypes['error']:
            msgstring[icode] = msgstring[icode] #+ '\n Execution continues.'
            if verbosity > 1: 
                pass
        else:
            if verbosity > 1: 
                pass
    elif icode == 0:
        for iicode in range(1,len(counter)):
            # printout the final stats for exception messages
            if not counter[iicode] == None:
                # print message as logger info 
                logger.warning(str(counter[iicode]) +' times '+ msgstring[iicode])
    return

exTypes = {'message': 'Message',
           'alert': 'Alert',
           'error': 'Error',
           'fatal': 'Fatal error'}

exceptions = [
    {'code': 1,
     'type': exTypes['fatal'],
     'msg': 'At least one required parameter is missing. Please enter all '
            'required inputs. Check documentation for details.'},
    {'code': 2,
     'type': exTypes['error'],
     'msg': 'zm (measurement height) must be larger than zero.'},
    {'code': 3,
     'type': exTypes['error'],
     'msg': 'z0 (roughness length) must be larger than zero.'},
    {'code': 4,
     'type': exTypes['error'],
     'msg': 'h (ABL height) must be larger than 10 m.'},
    {'code': 5,
     'type': exTypes['error'],
     'msg': 'zm (measurement height) must be smaller than h (PBL height).'},
    {'code': 6,
     'type': exTypes['alert'],
     'msg': 'zm (measurement height) should be above roughness sub-layer (12.5*z0).'},
    {'code': 7,
     'type': exTypes['error'],
     'msg': 'zm/ol (measurement height to Obukhov length ratio) must be equal or larger than -15.5.'},
    {'code': 8,
     'type': exTypes['error'],
     'msg': 'sigmav (standard deviation of crosswind) must be larger than zero.'},
    {'code': 9,
     'type': exTypes['error'],
     'msg': 'ustar (friction velocity) must be >=0.1.'},
    {'code': 10,
     'type': exTypes['error'],
     'msg': 'wind_dir (wind direction) must be >=0 and <=360.'},
    {'code': 11,
     'type': exTypes['fatal'],
     'msg': 'Passed data arrays (ustar, zm, h, ol) don\'t all have the same length.'},
    {'code': 12,
     'type': exTypes['fatal'],
     'msg': 'No valid zm (measurement height above displacement height) passed.'},
    {'code': 13,
     'type': exTypes['alert'],
     'msg': 'Using z0, ignoring umean if passed.'},
    {'code': 14,
     'type': exTypes['alert'],
     'msg': 'No valid z0 passed, using umean.'},
    {'code': 15,
     'type': exTypes['fatal'],
     'msg': 'No valid z0 or umean array passed.'},
    {'code': 16,
     'type': exTypes['error'],
     'msg': 'At least one required input is invalid. Skipping current footprint.'},
    {'code': 17,
     'type': exTypes['alert'],
     'msg': 'Only one value of zm passed. Using it for all footprints.'},
    {'code': 18,
     'type': exTypes['fatal'],
     'msg': 'if provided, rs must be in the form of a number or a list of numbers.'},
    {'code': 19,
     'type': exTypes['alert'],
     'msg': 'rs value(s) larger than 90% were found and eliminated.'},
    {'code': 20,
     'type': exTypes['error'],
     'msg': 'zm (measurement height) must be above roughness sub-layer (12.5*z0).'},
    {'code': 21,
     'type': exTypes['error'],
     'msg': 'umean (mean wind speed) must be >=0.0.'},
    ]

# =====================================================================================
# === plotting ===

def plot_footprint(cf):
    """
    Purpose:
    Usage:
    Side effects:
    Author: PRI/CE
    Date: Feb 2023
    """
    fpinfo = plot_footprint_parse_controlfile(cf)
    # read the netCDF file to a data structure
    ds = pfp_io.NetCDFRead(fpinfo["in_filename"])
    # check the required variables are in the data structure, if not may need to be calculated
    if not plot_footprint_check_variables(ds, fpinfo):
        return
    # update the info dictionary
    fpinfo["start"] = ds.root["Attributes"]["time_coverage_start"]
    fpinfo["end"] = ds.root["Attributes"]["time_coverage_end"]
    fpinfo["site_name"] = ds.root["Attributes"]["site_name"]
    # plot all the data
    plot_footprint_all(ds, fpinfo)
    # plot the seasonal footprint
    plot_footprint_seasonal(ds, fpinfo)
    return

def plot_footprint_all(ds, fpinfo):
    """
    Purpose:
    Usage:
    Side effects:
    Author: PRI/CE
    Date: Feb 2023
    """
    nbins = fpinfo["nbins"]
    wbins = fpinfo["wbins"]
    nsectors = fpinfo["nsectors"]
    site_name = ds.root["Attributes"]["site_name"]
    level = ds.root["Attributes"]["processing_level"]
    title = site_name + ": " + fpinfo["start"] + " to " + fpinfo["end"]
    title += "; all data"
    # read the variables
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Fsd = pfp_utils.GetVariable(ds, fpinfo["Fsd"])
    mask = numpy.ma.getmaskarray(Ws["Data"])
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Wd["Data"]))
    Ws["All"] = numpy.ma.masked_where(mask == True, Ws["Data"])
    Wd["All"] = numpy.ma.masked_where(mask == True, Wd["Data"])
    Ws["All"] = numpy.ma.compressed(Ws["All"])
    Wd["All"] = numpy.ma.compressed(Wd["All"])
    # prepare day/night data
    mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(Fsd["Data"]))
    for item in [Ws, Wd, Fsd]:
        item["DayNight"] = numpy.ma.masked_where(mask == True, item["Data"])
        item["DayNight"] = numpy.ma.compressed(item["DayNight"])
    ones = numpy.ones(len(Fsd["DayNight"]))
    zeros = numpy.zeros(len(Fsd["DayNight"]))
    idx = numpy.where(Fsd["DayNight"] <= fpinfo["Fsd_threshold"], ones, zeros)
    plt.ion()
    fig = plt.figure(figsize=(8, 8))
    fig.canvas.manager.set_window_title("footprint for "+fpinfo["site_name"]+", all data")
    gs = gridspec.GridSpec(2, 2) # 2 rows, 2 columns
    ax1=fig.add_subplot(gs[0, :], projection="footprint") # First row, first column
    ax2=fig.add_subplot(gs[1, 0], projection="footprint") # Second row, first column
    ax3=fig.add_subplot(gs[1, 1], projection="footprint") # Second row, second column
    ax1.box(Wd["All"], Ws["All"], bins=numpy.arange(0, nbins, wbins), nsector=nsectors)
    ax1.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25), fontsize=6)
    ax1.set_title("All")
    ax2.box(Wd["DayNight"][idx==0], Ws["DayNight"][idx==0],
            bins=numpy.arange(0, nbins, wbins),
            nsector=nsectors)
    ax2.set_title("Day")
    ax2.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25), fontsize=6)
    ax3.box(Wd["DayNight"][idx==1], Ws["DayNight"][idx==1],
            bins=numpy.arange(0, nbins, wbins),
            nsector=nsectors)
    ax3.set_title("Night")
    ax3.legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25), fontsize=6)
    fig.suptitle(title)
    fig.tight_layout()
    plt.draw()
    file_name = site_name.replace(" ","") + "_" + level + "_" + "footprint_all" + ".png"
    file_name = os.path.join(fpinfo["plot_path"], file_name)
    fig.savefig(file_name, format="png")
    pfp_utils.mypause(0.5)
    plt.ioff()
    return

def plot_footprint_seasonal(ds, fpinfo):
    """
    Purpose:
     Plot footprints for all day time and all night time data.
    Usage:
    Side effects:
    Author: PRI/CE
    Date: Feb 2023
    """
    nbins = fpinfo["nbins"]
    nsectors = fpinfo["nsectors"]
    wbins = fpinfo["wbins"]
    site_name = ds.root["Attributes"]["site_name"]
    level = ds.root["Attributes"]["processing_level"]
    title = site_name + ": " + fpinfo["start"] + " to " + fpinfo["end"]
    title += "; Seasonal"
    # read the variables
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    ldt["Month"] = numpy.array([dt.month for dt in ldt["Data"]])
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    seasons = {"Summer": [12, 1, 2], "Autumn": [3, 4, 5],
               "Winter": [6, 7, 8], "Spring": [9, 10, 11]}
    axs = {}
    plt.ion()
    fig = plt.figure(figsize=(8, 8))
    fig.canvas.manager.set_window_title("footprint for "+fpinfo["site_name"]+", seasons")
    for s, season in enumerate(seasons):
        cind = numpy.zeros(len(Ws["Data"]))
        for month in seasons[season]:
            idx = numpy.where(ldt["Month"] == month)[0]
            cind[idx] = 1
        wd = Wd["Data"][cind == 1]
        ws = Ws["Data"][cind == 1]
        mask = numpy.ma.getmaskarray(wd)
        mask = numpy.ma.mask_or(mask, numpy.ma.getmaskarray(ws))
        wd = numpy.ma.masked_where(mask == True, wd)
        wd = numpy.ma.compressed(wd)
        ws = numpy.ma.masked_where(mask == True, ws)
        ws = numpy.ma.compressed(ws)
        axs[season] = fig.add_subplot(2, 2, s+1, projection="footprint")
        axs[season].box(wd, ws,
                        bins=numpy.arange(0, nbins, wbins),
                        nsector=nsectors)
        axs[season].set_title(season)
        axs[season].legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25), fontsize=6)
    fig.suptitle(title)
    fig.tight_layout()
    plt.draw()
    file_name = site_name.replace(" ","") + "_" + level + "_" + "footprint_seasonal" + ".png"
    file_name = os.path.join(fpinfo["plot_path"], file_name)
    fig.savefig(file_name, format="png")
    pfp_utils.mypause(0.5)
    plt.ioff()
    return

# === utilities ===

def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    progress = round(progress,2)
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def crosswind_std(ds, fpinfo):
    msg = "No cross wind standard deviation in data set, Vsd will be calcuated from wind speed"
    logger.info(msg)
    logger.info(' Calculating cross wind standard deviation')
    # create a variable dictionary for Vsd
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    Vsd = pfp_utils.CreateEmptyVariable("Vsd", nrecs, datetime=ldt["Data"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    # get the required wind speed variables
    Ws = pfp_utils.GetVariable(ds, "Ws")
    # get the multiplier for calculating cross wind standard deviation from measured wind speed
    try:
        mult = fpinfo["Vsd_multiplier"]
    except:
        mult = 0.5   # default value is 50 % of measured wind speed
    Vsd["Data"] = mult * Ws["Data"]
    # get the QC flag
    Vsd["Flag"] = numpy.where(numpy.ma.getmaskarray(Vsd["Data"]) == True, ones, zeros)
    # update the variable attributes
    Vsd["Attr"]["units"] = "m/s"
    Vsd["Attr"]["long_name"] = "cross wind standard deviation"
    # put the cross wind standard deviation variable in the data structure
    pfp_utils.CreateVariable(ds, Vsd)
    return

def z0calc(ds, fpinfo):
    # aerodynamic roughness length, psi functions according to Dyer (1974)
    logger.info(' Calculating aerodynamic roughness length, Kormann and Meixner (2001) (Eqs. 31 to 35)')
    # get the required variables
    LM = pfp_utils.GetVariable(ds, fpinfo["ol"])
    Ws = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])
    # create a variable dictionary for z0
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    z0 = pfp_utils.CreateEmptyVariable("z0", nrecs, datetime=ldt["Data"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    #     a) create positive and negative LM masks
    LMp = numpy.ma.masked_where(LM["Data"] <  float(0),LM["Data"])
    LMn = numpy.ma.masked_where(LM["Data"] >= float(0),LM["Data"])
    # Calculate z0 assuming logarithmic wind profile, functions from Kormann and Meixner (2001) (Eqs. 31 to 35)
    #     b) for stable conditions, linear
    FIp = 5.0 * fpinfo["zm"]/LMp
    #     c) for unstable conditions
    zeta = (1.0-16.0*fpinfo["zm"]/LMn)**(0.25)
    FIn = -2.0*numpy.log(0.5*(1.0+zeta))-numpy.log(0.5*(1.0+zeta*zeta))+2.0*numpy.arctan(zeta)-0.5*c.Pi
    #     d) put both parts together again
    #FI = numpy.ma.mask_or(FIp,FIn)
    #     d1) fill positive and negative Fn masks
    FIp = numpy.ma.filled(FIp,float(0))
    FIn = numpy.ma.filled(FIn,float(0))
    FI  = FIp+FIn
    #     e) determine
    alpha = Ws["Data"] * 0.4 / ustar["Data"] - FI
    #     f) finally calculate the roughness length
    z0["Data"] = fpinfo["zm"] / numpy.exp(alpha)
    #set a lower limit for z0 to avoid numeric problems
    z0["Data"] = numpy.ma.masked_where(z0["Data"]<0.0001,z0["Data"])
    z0["Data"] = numpy.ma.filled(z0["Data"],0.0001)
    # get the QC flag
    z0["Flag"] = numpy.where(numpy.ma.getmaskarray(z0["Data"]) == True, ones, zeros)
    # update the variable attributes
    z0["Data"]["Attr"]["units"] = "m"
    z0["Data"]["Attr"]["long_name"] = "roughness lenth"
    # put the cross wind standard deviation variable in the data structure
    pfp_utils.CreateVariable(ds, z0)
    return
