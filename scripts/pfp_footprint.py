import logging
import os
import copy
import warnings
# 3rd party
import matplotlib
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    matplotlib.rcParams['toolbar'] = 'toolmanager'
import numpy
# PFP modules
from scripts import constants as c
from scripts import pfp_classes
from scripts import pfp_ck
from scripts import pfp_io
from scripts import pfp_utils
from scripts import pfp_ts

"""
Derive a flux footprint estimate based on the paper from Kormann and Meixner (2001). The equations solved
follow the work from Neftel et al. (2008) and the ART footprint tool. Currently no sections are 
defined to estimate the proportion of footprint source from this defined areas as in Neftel et al. (2008).
See Kormann, R. and Meixner, F.X., 2001: An analytical footprint model for non-neutral stratification.
Boundary-Layer Meteorology 99: 207. https://doi.org/10.1023/A:1018991015119 and 
Neftel, A., Spirig, C., Ammann, C., 2008: Application and test of a simple tool for operational footprint 
evaluations. Environmental Pollution 152, 644-652. for details.

FKM Input
    zm       = Measurement height above displacement height (i.e. z-d) [m]
    z0       = Roughness length [m], calculated according to Kormann and Meixner, 2001
    umean    = Vector of mean wind speed at zm [ms-1]
    ol       = Vector of Obukhov length [m]
    Vsd      = Vector of standard deviation of lateral velocity fluctuations [ms-1]
               if not available it is calculated as 0.5*Ws (wind speed)
    ustar    = Vector of friction velocity [ms-1]
    Ws       = Vector of wind direction in degrees (of 360) for rotation of the footprint

Footprint calculation following the Kormann-Meixner Approach
1st Version by Marx Stampfli, stampfli mathematics, Bern, Switzerland
Revision December 2006, C. Spirig and A. Neftel, Agroscope, Zurich, Switzerland
            ported into python CM Ewenz, Adelaide Nov 2015 to May 2016
Array of lower and upper ranges for following parameters
        ustar,     LM,  Vsd,  Wd,   zm,    Ws
lrange = 0.01,-999999,    0,   0,    0,     0
urange =    5, 999999,   20, 360, 1000,    30

Created: May 2018 Cacilia Ewenz
         version: 0.1
Changes: 08/06/2018 Cacilia Ewenz
         22/02/2023 Calculate Kormann-Meixner parameters: A, B, C, D, E, x0, KM_p01a, KM_p01b, KM_p01c
         
Copyright (C) 2018, 2019, 2020, 2021, 2022, 2023, Cacilia Ewenz
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
    in_filename = pfp_io.get_infilenamefromcf(cf)
    if not pfp_utils.file_exists(in_filename):
        logger.error(" Unable to find netCDF file " + in_filename)
        return
    ds = pfp_io.NetCDFRead(in_filename)
    if ds.info["returncodes"]["value"] != 0:
        return
    # read the controlfile
    fpinfo = calculate_footprint_parse_controlfile(cf)
    # read the netCDF file to a data structure
    ds_fkm = pfp_classes.DataStructure()
    # copy global attributes to new dictionary
    ds_fkm = copy.deepcopy(ds)
    # check the required variables are in the data structure
    if not calculate_footprint_check_variables(ds_fkm, fpinfo):
        return
    # (1) Evaluate input data and apply constraints
    logger.info("Check footprint input variables")    
    check_fkm_inputs(ds_fkm,fpinfo)
    # (2) get a mask for all input variables
    logger.info("Mask footprint input variables to valids")    
    composite_mask(ds_fkm, fpinfo)
    # (3) calculate Kormann and Meixner, 2001 footprint
    logger.info("Calculate Kormann-Meixner footprint variables")    
    kormann_meixner(ds_fkm, fpinfo)
    # (4) Do range checks on output parameter
    # apply the quality control checks (range, diurnal, exclude dates and exclude hours
    #pfp_ck.do_qcchecks(cf, ds_fkm) or directly pfp_ck.do_rangecheck(cf,ds,section,series,code=2)
    # (5) write data to level 3 output netCDF file
    logger.info("Write output to L3 footprint file")        
    out_filename = fpinfo["out_filename"]
    pfp_io.NetCDFWrite(out_filename, ds_fkm)
    logger.info("Finished Footprint calculation")
    
    return

def calculate_footprint_check_variables(ds, fpinfo):
    logger.info("Check the variables named in the control file exist in the data structure.")
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
    # get the input and output file names from the control file
    fpinfo["in_filename"] = pfp_io.get_infilenamefromcf(cf)
    try:
        # everything else gets here
        path = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "file_path")
        name = pfp_utils.get_keyvaluefromcf(cf, ["Files"], "out_filename")
        fpinfo["out_filename"] = os.path.join(path, name)
    except:
        fpinfo["out_filename"] = fpinfo["in_filename"].replace(".nc", "_fp.nc")
    # get the day/night separator Fsd name and its threshold
    opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "Fsd_threshold", default="10")
    fpinfo["Fsd_threshold"] = int(opt)
    # get the wind direction and wind speed variable names
    fpinfo["Wd"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Wd"], "name", default="Wd")
    fpinfo["Ws"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Ws"], "name", default="Ws")
    fpinfo["Fsd"] = pfp_utils.get_keyvaluefromcf(cf, ["Variables", "Fsd"], "name", default="Fsd")
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
    ol    = pfp_utils.GetVariable(ds, fpinfo["ol"])
    Ws    = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd    = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd   = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])
    
    # get the instrument height above the displacement height (2/3*canopy_height)
    th = float(pfp_utils.strip_non_numeric(ds.root["Attributes"]["tower_height"]))
    dh = float(pfp_utils.strip_non_numeric(ds.root["Attributes"]["canopy_height"]))
    fpinfo["zm"] = th - (2/3)*dh
    # get the roughness length and z/L
    z0calc(ds, fpinfo)
    z0 = pfp_utils.GetVariable(ds, "z0")
    
    # create a variable dictionary for zmol
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    # Create zm/L variable
    zmol = pfp_utils.CreateEmptyVariable("zmol", nrecs, datetime=ldt["Data"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    zmol["Data"] = fpinfo["zm"]/ol["Data"]
    # get the QC flag
    zmol["Flag"] = numpy.where(numpy.ma.getmaskarray(zmol["Data"]) == True, ones, zeros)
    # update the variable attributes
    zmol["Attr"]["units"] = "1"
    zmol["Attr"]["long_name"] = "zm over L"
    # put the cross wind standard deviation variable in the data structure
    pfp_utils.CreateVariable(ds, zmol)
    
    # start checking constraints for footprint calculation
    if fpinfo["zm"] <= 0.:
        msg = "zm (measurement height) must be larger than zero."
        logger.error(msg)
        return False
    z0["Data"]    = numpy.ma.masked_less(z0["Data"], 0.0, copy=True)
    zmol["Data"]  = numpy.ma.masked_less(zmol["Data"], -3.0, copy=True)
    zmol["Data"]  = numpy.ma.masked_greater(zmol["Data"], 3.0, copy=True)
    ol["Data"]    = numpy.ma.masked_less(ol["Data"], -99999.0, copy=True)
    ol["Data"]    = numpy.ma.masked_greater(ol["Data"], 99999.0, copy=True)
    Vsd["Data"]   = numpy.ma.masked_less_equal(Vsd["Data"], 0.0, copy=True)
    Vsd["Data"]   = numpy.ma.masked_greater_equal(Vsd["Data"], 20.0, copy=True)
    ustar["Data"] = numpy.ma.masked_less_equal(ustar["Data"], 0.1, copy=True)
    ustar["Data"] = numpy.ma.masked_greater_equal(ustar["Data"], 5.0, copy=True)
    Ws["Data"]    = numpy.ma.masked_less_equal(Ws["Data"], 0.1, copy=True)
    Ws["Data"]    = numpy.ma.masked_greater_equal(Ws["Data"], 30.0, copy=True)
    Wd["Data"]    = numpy.ma.masked_greater(Wd["Data"], 360.0, copy=True)
    Wd["Data"]    = numpy.ma.masked_less(Wd["Data"], 0.0, copy=True)
    return True

def composite_mask(ds, fpinfo):
    """
    Purpose:
     Calculate and apply a composite mask for all Kormann-Maixner input variables.
    Author: CE
    Date: Feb 2023
    """
    ldt   = pfp_utils.GetVariable(ds, "DateTime")
    # get the required variables
    z0    = pfp_utils.GetVariable(ds, "z0")
    zmol  = pfp_utils.GetVariable(ds, "zmol")
    ol    = pfp_utils.GetVariable(ds, fpinfo["ol"])
    Ws    = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd    = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd   = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])

    # get a composite mask over all variables needed for Kormann-Meixner 2001 footprint
    mask_all = numpy.ma.getmaskarray(ds.root["Variables"]["ustar"])
    for item in ["Ws", "ol", "Vsd", "Wd", "z0", "zmol"]:
        mask_item = numpy.ma.getmaskarray(item)
        mask_all = numpy.ma.mask_or(mask_all, mask_item)
    
    # apply composite mask but don't remove masked elements
    Ws["Data"]    = numpy.ma.masked_where(mask_all == True, Ws["Data"])
    ol["Data"]    = numpy.ma.masked_where(mask_all == True, ol["Data"])
    zmol["Data"]  = numpy.ma.masked_where(mask_all == True, zmol["Data"])
    Vsd["Data"]   = numpy.ma.masked_where(mask_all == True, Vsd["Data"])
    ustar["Data"] = numpy.ma.masked_where(mask_all == True, ustar["Data"])
    Wd["Data"]    = numpy.ma.masked_where(mask_all == True, Wd["Data"])
    z0["Data"]    = numpy.ma.masked_where(mask_all == True, z0["Data"])
    ldt["Data"]   = numpy.ma.masked_where(mask_all == True, ldt["Data"])
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
    ldt   = pfp_utils.GetVariable(ds, "DateTime")
    nrecs = len(ldt["Data"])
    # get the required variables
    z0    = pfp_utils.GetVariable(ds, "z0")
    ol    = pfp_utils.GetVariable(ds, fpinfo["ol"])
    zmol  = pfp_utils.GetVariable(ds, "zmol")
    Ws    = pfp_utils.GetVariable(ds, fpinfo["Ws"])
    Wd    = pfp_utils.GetVariable(ds, fpinfo["Wd"])
    Vsd   = pfp_utils.GetVariable(ds, fpinfo["Vsd"])
    ustar = pfp_utils.GetVariable(ds, fpinfo["ustar"])
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    # calculate height above displacement height
    zm = fpinfo["zm"]
    # a) create positive, zero and negative LM masks
    LMp = numpy.ma.masked_where(ol["Data"] <  float(0),ol["Data"])
    LMz = numpy.ma.masked_where(ol["Data"] == float(0),ol["Data"])
    LMn = numpy.ma.masked_where(ol["Data"] > float(0),ol["Data"])
    # phi, u, m, n
    p_phi = 1.0+5.0*zmol["Data"]
    p_u   = (ustar["Data"] / c.k) * (numpy.log(zm / z0["Data"]) + 5.0 * zmol["Data"])
    p_m   = (1.0 + 5.0 * (zmol["Data"])) / (numpy.log(zm / z0["Data"]) + 5.0 * zmol["Data"])
    p_n   = 1.0 / (1.0 + 5.0 * zmol["Data"])
    z_phi = ones
    z_u   = ustar["Data"] / c.k * (numpy.log(zm / z0["Data"]))
    z_m   = ustar["Data"] / c.k / z_u
    z_n   = ones
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
    Kmaj = c.k * ustar["Data"] * zm / phi / zm**n
    # Kmaj corresponds to kappa in KM 2001
    xi = Umaj * (zm **r) / (r * r * Kmaj)
    # xPhiMax is the (x-)position of the maximum of phi
    xPhiMax = r * xi / (2 * r + 1)
    # determine GammaProxmu and GammaProx1r depending on mu and 1/r
    GammaProxmu = zeros
    GammaProxmu = numpy.ma.masked_where(mu == 0.0, GammaProxmu)
    GammaProxmu = (1.0 / mu) + 0.1002588 * numpy.exp(mu) - 0.493536 + 0.3066 * mu - 0.09 * (mu**2)
    GammaProxmu  = numpy.ma.filled(GammaProxmu,float(0))
    GammaProx1r = zeros
    GammaProx1r = numpy.ma.masked_where(1/r == 0.0, GammaProx1r)
    GammaProx1r = r + 0.1002588 * numpy.exp(1/r) - 0.493536 + 0.3066 * (1/r) - 0.09 * ((1/r)**2)
    GammaProx1r = numpy.ma.filled(GammaProx1r,float(0))
    # Kormann-Meixner parameters A-E
    A = pfp_utils.CreateEmptyVariable("A", nrecs, datetime=ldt["Data"])
    B = pfp_utils.CreateEmptyVariable("B", nrecs, datetime=ldt["Data"])
    C = pfp_utils.CreateEmptyVariable("C", nrecs, datetime=ldt["Data"])
    D = pfp_utils.CreateEmptyVariable("D", nrecs, datetime=ldt["Data"])
    E = pfp_utils.CreateEmptyVariable("E", nrecs, datetime=ldt["Data"])
    KM_p01a = pfp_utils.CreateEmptyVariable("KM_p01a", nrecs, datetime=ldt["Data"])
    KM_p01b = pfp_utils.CreateEmptyVariable("KM_p01b", nrecs, datetime=ldt["Data"])
    #KM_p01c = pfp_utils.CreateEmptyVariable("KM_p01c", nrecs, datetime=ldt["Data"])
    A["Data"] = 1 + mu
    B["Data"] = Umaj * (zm**r) / r / r / Kmaj
    C["Data"] = (B["Data"]**mu) /GammaProxmu
    D["Data"] = Vsd["Data"] * GammaProx1r / GammaProxmu / ((r * r * Kmaj / Umaj)**(m / r)) / Umaj
    E["Data"] = (r - m) / r
    # Ellipse parameters: Calculate the near and far distance for 1% of the max level
    philevel = 0.01
    x0 = xPhiMax / 3
    p = philevel * xPhiMax
    for _ in [0,1,2,3,4,5]:
        Fm = C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"]) - p
        dFm = (C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"])) / x0**2 * -(A["Data"] * x0 - B["Data"]) - p
        x0 = x0 - Fm / dFm
    KM_p01a["Data"] = x0
    x0 = xPhiMax * 3
    p = philevel * xPhiMax
    for _ in [0,1,2,3,4,5]:
        Fm = C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"]) - p
        dFm = (C["Data"] * numpy.exp(-B["Data"] / x0) * x0**(-A["Data"])) / x0**2 * -(A["Data"] * x0 - B["Data"]) - p
        x0 = x0 - Fm / dFm
    KM_p01b["Data"] = x0
    # set x0 = center of ellipse
    #x0 = 0.5*(KM_p01b["Data"] - KM_p01a["Data"])
    # Delivers y for which phi(xcenter,y)=level, i.e. "width of ellipse" ellipse parameter d: short axes of ellipse
    #for _ in [0,1,2,3,4,5]:
    #    help1 = D["Data"] * x0**E["Data"] * numpy.sqrt(2)
    #    help2 = C["Data"] * numpy.exp(-B["Data"] / x0)
    #    help3 = x0**(A["Data"] + E["Data"])
    #    y = help1 * numpy.sqrt(numpy.log(help2 / numpy.sqrt(2 * numpy.pi) / help3 / D["Data"] / p))
    #    help4 = (1 - y**2 / D["Data"]**2 / x0**(2 * E["Data"]))
    #    Gm = help4 * E["Data"] / x0 - xi / x0**2 + (1 + mu) / x0
    #    help5 = (y**2 * E["Data"]**2)
    #    help6 = x0**(2 * E["Data"] + 2)
    #    help7 = (1 - y**2 / D["Data"]**2 / x0**(2 * E["Data"]))
    #    dGm = 2 * help5 / D["Data"]**2 / help6 - help7 * E["Data"] / x0**2 + 2 * xi / x0**3 + (1 + mu) / x0**2
    #    x0 = x0 - Gm / dGm
       
    #KM_p01c["Data"] = D["Data"] * x0**E["Data"] * numpy.sqrt(2) * numpy.sqrt(numpy.log(C["Data"] * numpy.exp(-B["Data"] / x0) / numpy.sqrt(2 * numpy.pi) / x0**(A["Data"] + E["Data"]) / D["Data"] / p))
    # range check for ellipse parameter
    KM_p01a["Data"]    = numpy.ma.masked_greater(KM_p01a["Data"],  10000.0, copy=True)
    KM_p01b["Data"]    = numpy.ma.masked_greater(KM_p01b["Data"], 100000.0, copy=True)
    #KM_p01c["Data"]    = numpy.ma.masked_greater(KM_p01c["Data"],  50000.0, copy=True)
    # get the QC flag
    A["Flag"] = numpy.where(numpy.ma.getmaskarray(A["Data"]) == True, ones, zeros)
    B["Flag"] = numpy.where(numpy.ma.getmaskarray(B["Data"]) == True, ones, zeros)
    C["Flag"] = numpy.where(numpy.ma.getmaskarray(C["Data"]) == True, ones, zeros)
    D["Flag"] = numpy.where(numpy.ma.getmaskarray(D["Data"]) == True, ones, zeros)
    E["Flag"] = numpy.where(numpy.ma.getmaskarray(E["Data"]) == True, ones, zeros)
    KM_p01a["Flag"] = numpy.where(numpy.ma.getmaskarray(KM_p01a["Data"]) == True, ones, zeros)
    KM_p01b["Flag"] = numpy.where(numpy.ma.getmaskarray(KM_p01b["Data"]) == True, ones, zeros)
    #KM_p01c["Flag"] = numpy.where(numpy.ma.getmaskarray(KM_p01c["Data"]) == True, ones, zeros)
    # update the variable attributes
    A["Attr"]["units"] = "1"
    B["Attr"]["units"] = "1"
    C["Attr"]["units"] = "1"
    D["Attr"]["units"] = "1"
    E["Attr"]["units"] = "1"
    KM_p01a["Attr"]["units"] = "m"
    KM_p01b["Attr"]["units"] = "m"
    #KM_p01c["Attr"]["units"] = "m"
    A["Attr"]["long_name"] = "Kormann-Meixner parameter A"
    B["Attr"]["long_name"] = "Kormann-Meixner parameter B"
    C["Attr"]["long_name"] = "Kormann-Meixner parameter C"
    D["Attr"]["long_name"] = "Kormann-Meixner parameter D"
    E["Attr"]["long_name"] = "Kormann-Meixner parameter E"
    KM_p01a["Attr"]["long_name"] = "Kormann-Meixner intersection point of 1%- ellipse nearest to sensor"
    KM_p01b["Attr"]["long_name"] = "Kormann-Meixner intersection point of 1%- ellipse farthest from sensor"
    #KM_p01c["Attr"]["long_name"] = "Kormann-Meixner half width of 1%-ellipse"
    pfp_utils.CreateVariable(ds, A)
    pfp_utils.CreateVariable(ds, B)
    pfp_utils.CreateVariable(ds, C)
    pfp_utils.CreateVariable(ds, D)
    pfp_utils.CreateVariable(ds, E)
    pfp_utils.CreateVariable(ds, KM_p01a)
    pfp_utils.CreateVariable(ds, KM_p01b)
    #pfp_utils.CreateVariable(ds, KM_p01c)
    return

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
    logger.info('Calculating aerodynamic roughness length, Kormann and Meixner (2001) (Eqs. 31 to 35)')
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
    z0["Attr"]["units"] = "m"
    z0["Attr"]["long_name"] = "roughness lenth"
    # put the cross wind standard deviation variable in the data structure
    pfp_utils.CreateVariable(ds, z0)
    return
