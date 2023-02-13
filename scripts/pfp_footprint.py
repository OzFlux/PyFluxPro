import datetime
import logging
import math
import os
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
from scripts import pfp_ck
from scripts import pfp_io
from scripts import pfp_utils
from scripts import meteorologicalfunctions as pfp_mf

logger = logging.getLogger("pfp_log")

def calculate_footprint(cf):
    """
    Purpose:
    Usage:
    Side effects:
    Author: PRI/CE
    Date: Feb 2023
    """
    fpinfo = calculate_footprint_parse_controlfile(cf)
    # read the netCDF file to a data structure
    ds = pfp_io.NetCDFRead(fpinfo["in_filename"])
    # check the required variables are in the data structure
    if not calculate_footprint_check_variables(ds, fpinfo):
        msg = 
    # update the info dictionary
    fpinfo["start"] = ds.root["Attributes"]["time_coverage_start"]
    fpinfo["end"] = ds.root["Attributes"]["time_coverage_end"]
    fpinfo["site_name"] = ds.root["Attributes"]["site_name"]
    # plot all the data
    plot_footprint_all(ds, fpinfo)
    # plot the seasonal footprint
    plot_footprint_seasonal(ds, fpinfo)
    return

def calculate_footprint_check_variables(ds, fpinfo):
    """ Check the variables named in the control file exist in the data structure."""
    labels = list(ds.root["Variables"].keys())
    for item in ["Fsd", "Wd", "Ws", "ustar", "Vsd", "ol"]:
        name = fpinfo[item]
        if name not in labels:
            msg = name + "(" + item + ") not found in " + fpinfo["in_filename"]
            logger.error(msg)
            ok = False
            if name in ["ustar", "Vsd", "ol"]:
                msg = name + "(" + item + ") will be calculated."
                logger.error(msg)
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
    # get the input file name from the control file
    #fpinfo["in_filename"] = pfp_io.get_infilenamefromcf(cf)
    #opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "number_speed_bins", default="8")
    #fpinfo["nbins"] = int(opt)
    #opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "speed_bin_width", default="1")
    #fpinfo["wbins"] = int(opt)
    #opt = pfp_utils.get_keyvaluefromcf(cf, ["Options"], "number_sectors", default="16")
    #fpinfo["nsectors"] = int(opt)
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

def CalculateMoninObukhovLength(ds, fpinfo):
    """
    Purpose:
     Calculate the Monin Obukhov length.
    Usage:
     CalculateMoninObukhovLength(ds)
     where ds is a data structure
    Side effects:
     Creates a new series in the data structure containing the Monin-Obukhov length.
    Author: PRI
    Date: April 2018
    """
    logger.info(' Calculating Monin-Obukhov length')
    # create a variable dictionary for L
    nrecs = int(ds.root["Attributes"]["nc_nrecs"])
    ldt = pfp_utils.GetVariable(ds, "DateTime")
    L = pfp_utils.CreateEmptyVariable("L", nrecs)
    # create QC flags
    zeros = numpy.zeros(nrecs, dtype=numpy.int32)
    ones = numpy.ones(nrecs, dtype=numpy.int32)
    # get the required meteorological variables
    Ta = pfp_utils.GetVariable(ds, "Ta")
    ps = pfp_utils.GetVariable(ds, "ps")
    if "e" in list(ds.series.keys()): 
        vp = pfp_utils.GetVariable(ds, "e")
    else:
        vp = pfp_utils.CreateEmptyVariable("vp", nrecs)
        Ah = pfp_utils.GetVariable(ds, "Ah")
        vp["Data"] = pfp_mf.vapourpressure(Ah["Data"],Ta["Data"])
    
    # get the required fluxes
    ustar = pfp_utils.GetVariable(ds, "ustar")
    Fh    = pfp_utils.GetVariable(ds, "Fh")
    # calculate the density of dry air
    rho_dry = pfp_mf.densitydryair(Ta["Data"], ps["Data"], vp["Data"])
    # calculate virtual potential temperature
    Tp = pfp_mf.theta(Ta["Data"], ps["Data"])
    mr = pfp_mf.mixingratio(ps["Data"], vp["Data"])
    Tvp = pfp_mf.virtualtheta(Tp, mr)
    L["Data"] = -Tvp*rho_dry*c.Cp*(ustar["Data"]**3)/(c.g*c.k*Fh["Data"])
    # get the QC flag
    L["Flag"] = numpy.where(numpy.ma.getmaskarray(L["Data"]) == True, ones, zeros)
    # update the variable attributes
    L["Attr"]["units"] = "m"
    L["Attr"]["long_name"] = "Monin-Obukhov length"
    L["Attr"]["standard_name"] = "not defined"
    # put the Monin-Obukhov variable in the data structure
    pfp_utils.CreateVariable(ds, L)
    return

def crosswind_std(ds, fpinfo):
    msg = "No cross wind standard deviation in data set, Vsd will be calcuated from wind speed"
    logger.info(msg)
    
    return

def z0calc(ds, fpinfo):
    # aerodynamic roughness length
    # Psi functions according to Dyer (1974)
    #     a) create positive and negative LM masks
    LMp = numpy.ma.masked_where(LM <  float(0),LM)
    LMn = numpy.ma.masked_where(LM >= float(0),LM)
    # Calculate z0 assuming logarithmic wind profile
    #          === functions are from Kormann and Meixner (2001) (Eqs. 31 to 35)
    #     b) for stable conditions, linear
    FIp = 5.0 * zm/LMp
    #     c) for unstable conditions
    zeta = (1.0-16.0*zm/LMn)**(0.25)
    FIn = -2.0*numpy.log(0.5*(1.0+zeta))-numpy.log(0.5*(1.0+zeta*zeta))+2.0*numpy.arctan(zeta)-0.5*c.Pi
    #     d) put both parts together again
    #FI = numpy.ma.mask_or(FIp,FIn)
    #     d1) fill positive and negative Fn masks
    FIp = numpy.ma.filled(FIp,float(0))
    FIn = numpy.ma.filled(FIn,float(0))
    FI  = FIp+FIn
    #     e) determine
    alpha = U_meas * 0.4 / UStar - FI
    #     f) finally calculate the roughness length
    ZNull = zm / numpy.exp(alpha)
    #!#            === functions derived from Leclerc and Foken 2015 book, page 61 after Hogstroem, 1988
    #!#     b) for stable conditions, linear
    #!FIp = -6.0 * zm/LMp
    #!#     c) for unstable conditions
    #!zeta = (1.0+19.3*zm/LMn)**0.25
    #!temp = 0.125*(1.0+zeta*zeta)*(1.0+zeta)*(1.0+zeta)
    #!FIn = numpy.log(temp)-2.0*numpy.arctan(zeta)+0.5*c.Pi
    #!#     d) put both parts together again
    #!#FI = numpy.ma.mask_or(FIp,FIn,copy=True)
    #!# d1) fill positive and negative Fn masks
    #!FIp = numpy.ma.filled(FIp,float(0))
    #!FIn = numpy.ma.filled(FIn,float(0))
    #!FI  = FIp+FIn
    #!#     e) determine
    #!alpha = U_meas * 0.4 / UStar + FI
    #!#     f) finally calculate the roughness length
    #!ZNull = zm / numpy.exp(alpha)
    #!#            ===
    #set a lower limit for z0 to avoid numeric problems
    ZNull = numpy.ma.masked_where(ZNull<0.0001,ZNull)
    ZNull = numpy.ma.filled(ZNull,0.0001)
    return ZNull

def BLH(ol, ustar, lat, zm):
    #def BLH(ol, ustar, dt, lat, zm):
    # --- if no boundary layer height available use Kljun et al., 2015 analytical solution for Habl
    #   blh for stable and neutral conditions - Nieuwstadt (1981)
    #            h = (L/3.8)*(-1 + (1 + 2.28*(ustar/(f*L)))^0.5)                          (Eq.1.1)
    #              with L = MOL, ustar - friction velocity, g = acceleration due to gravity and
    #                   f = 2 omega sin (phi) = Coriolis parameter
    #   blh for convective conditions (needs to be integrated due to near symmetric diurnal cycle
    #            of surface heat flux). The resulting rate of change for the boundary layer height
    #            is implicit in h and may be solved iteratively or using h(t_i) to determine the rate
    #            of change to yield h(t_i+1)
    #            at sunrise before Fh becomes positive - use Eq.1.1 for initial conditions
    #               then dh/dt = bar(w'Theta')_0 / gamma*[(h^2 / ((1+2*A)*h - 2*B*k*L))   (Eq.1.2)
    #                           + ((C * ustar^2 * T) / (gamma*g*[(1+A)*h - B*k*L]))]^(-1)
    #            gamma = gradient of potential temp above convective bl; ~0.01 K/m for typical midlatitude
    #                      (dh/dt quite sensitive to gamma)
    #            A = 0.2; B + 2.5; C = 8 (derived from similarity relations in cbl)
    #   so make sure H is availble from external sources, you really don't want to calculate it
    #A = 0.2
    #B = 2.5
    #C = 8
    #gamma = 0.01          # K/m
    omega = 0.000072921   # rad/s
    f = 2.0 * omega * numpy.sin(lat * numpy.pi / 180.)
    if ol > 0:    # stable
        Habl = (ol/3.8)*(-1 + 2.28*(ustar/(f*ol)))**0.5
    elif ol <= 0: # convective
        # here we need to model the data day by day
        Habl = (ol/3.8)*(-1 + 2.28*(ustar/(f*ol)))**0.5
    elif float(zm)/ol <= -15.5:
        # not valid in Kljun et al., 2015
        Habl = -9999
    return Habl

