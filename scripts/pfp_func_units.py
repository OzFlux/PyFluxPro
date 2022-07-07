# standard modules
import logging
# 3rd party
import numpy
# PFP modules
from scripts import meteorologicalfunctions as pfp_mf
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def fraction_to_percent(ds, RH_out, RH_in):
    """
    Purpose:
     Function to convert RH in units of "frac" (0 to 1) to "percent" (1 to 100).
    Usage:
     pfp_func_units.fraction_to_percent(ds, RH_out, RH_in)
    Author: PRI
    Date: August 2019
    """
    var_in = pfp_utils.GetVariable(ds, RH_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "percent", mode="quiet")
    var_out["Label"] = RH_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def gH2Opm3_to_percent(ds, RH_out, AH_in, Ta_in):
    """
    Purpose:
     Function to convert absolute humidity in units of g/m^3 to relative humidity in percent.
    Usage:
     pfp_func_units.gH2Opm3_to_percent(ds, RH_out, AH_in, Ta_in)
    Author: PRI
    Date: September 2020
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [AH_in, Ta_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + RH_out + " not calculated"
            logger.error(msg)
            return 0
    AH = pfp_utils.GetVariable(ds, AH_in)
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    RH = pfp_utils.GetVariable(ds, RH_out)
    RH["Data"] = pfp_mf.relativehumidityfromabsolutehumidity(AH["Data"], Ta["Data"])
    RH["Flag"] = numpy.where(numpy.ma.getmaskarray(RH["Data"]) == True, ones, zeros)
    RH["Attr"]["units"] = "percent"
    pfp_utils.CreateVariable(ds, RH)
    return 1

def gH2Opm3_to_mmolpm3(ds, H2O_out, AH_in):
    """
    Purpose:
     Calculate H2O molar density in mmol/m^3 from absolute humidity in g/m^3.
    Usage:
     pfp_func_units.gH2Opm3_to_mmolpm3(ds, MD_out, AH_in)
    Author: PRI
    Date: September 2020
    """
    for item in [AH_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + H2O_out + " not calculated"
            logger.error(msg)
            return 0
    var_in = pfp_utils.GetVariable(ds, AH_in)
    got_variance = False
    if var_in["Label"][-3:] == "_Vr" and var_in["Attr"]["units"] in ["g^2/m^6", "gH2O^2/m^6"]:
        got_variance = True
        var_in["Data"] = numpy.ma.sqrt(var_in["Data"])
        var_in["Attr"]["units"] = "g/m^3"
    var_out = pfp_utils.convert_units_func(ds, var_in, "mmol/m^3", mode="quiet")
    var_out["Label"] = H2O_out
    if got_variance:
        var_out["Data"] = var_out["Data"]*var_out["Data"]
        var_out["Attr"]["units"] = "mmol^2/m^6"
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def gH2Opm3_to_mmolpmol(ds, MF_out, AH_in, Ta_in, ps_in):
    """
    Purpose:
     Calculate H2O mole fraction in mml/mol from absolute humidity in g/m^3.
    Usage:
     pfp_func_units.gH2Opm3_to_mmolpmol(ds, MF_out, AH_in, Ta_in, ps_in)
    Author: PRI
    Date: August 2019
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [AH_in, Ta_in, ps_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + MF_out + " not calculated"
            logger.error(msg)
            return 0
    AH = pfp_utils.GetVariable(ds, AH_in)
    AH = pfp_utils.convert_units_func(ds, AH, "g/m^3")
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    Ta = pfp_utils.convert_units_func(ds, Ta, "degC")
    ps = pfp_utils.GetVariable(ds, ps_in)
    ps = pfp_utils.convert_units_func(ds, ps, "kPa")
    MF = pfp_utils.GetVariable(ds, MF_out)
    MF["Data"] = pfp_mf.h2o_mmolpmolfromgpm3(AH["Data"], Ta["Data"], ps["Data"])
    MF["Flag"] = numpy.where(numpy.ma.getmaskarray(MF["Data"]) == True, ones, zeros)
    MF["Attr"]["units"] = "mmol/mol"
    pfp_utils.CreateVariable(ds, MF)
    return 1

def hPa_to_kPa(ds, ps_out, ps_in):
    """
    Purpose:
     Function to convert pressure from hPa (mb) to kPa.
    Usage:
     pfp_func_units.hPa_to_kPa(ds, ps_in, ps_out)
    Author: PRI
    Date: February 2018
    """
    var_in = pfp_utils.GetVariable(ds, ps_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "kPa", mode="quiet")
    var_out["Label"] = ps_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def K_to_C(ds, T_out, T_in):
    """
    Purpose:
     Function to convert temperature from K to C.
    Usage:
     pfp_func_units.K_to_C(ds, T_out, T_in)
    Author: PRI
    Date: February 2018
    """
    if T_in not in list(ds.root["Variables"].keys()):
        msg = " Convert_K_to_C: variable " + T_in + " not found, skipping ..."
        logger.warning(msg)
        return 0
    if "<" in T_out or ">" in T_out:
        logger.warning(" ***")
        msg = " *** " + T_in + ": illegal name (" + T_out + ") in function, skipping ..."
        logger.warning(msg)
        logger.warning(" ***")
        return 0
    var_in = pfp_utils.GetVariable(ds, T_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "degC", mode="quiet")
    var_out["Label"] = T_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def kgpm3_to_gpm3(ds, AH_out, AH_in):
    """
    Purpose:
     Function to convert absolute humidity from kg/m^3 to g/m^3.
    Usage:
     pfp_func_units.kgpm3_to_gpm3(ds, AH_out, AH_in)
    Author: PRI
    Date: August 2020
    """
    var_in = pfp_utils.GetVariable(ds, AH_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "g/m^3", mode="quiet")
    var_out["Label"] = AH_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def mgCO2pm3_to_mmolpm3(ds, CO2_out, CO2_in):
    """
    Purpose:
     Calculate CO2 molar density in mmol/m3 from CO2 concentration in mg/m3.
    Usage:
     pfp_func_units.mgCO2pm3_to_mmolpm3(ds, CO2_out, CO2_in)
    Author: PRI
    Date: September 2020
    """
    for item in [CO2_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + CO2_out + " not calculated"
            logger.error(msg)
            return 0
    var_in = pfp_utils.GetVariable(ds, CO2_in)
    got_variance = False
    if var_in["Label"][-3:] == "_Vr" and var_in["Attr"]["units"] in ["mg^2/m^6", "mgCO2^2/m^6"]:
        got_variance = True
        var_in["Data"] = numpy.ma.sqrt(var_in["Data"])
        var_in["Attr"]["units"] = "mg/m^3"
    var_out = pfp_utils.convert_units_func(ds, var_in, "mmol/m^3", mode="quiet")
    var_out["Label"] = CO2_out
    if got_variance:
        var_out["Data"] = var_out["Data"]*var_out["Data"]
        var_out["Attr"]["units"] = "mmol^2/m^6"
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def mgCO2pm3_to_umolpmol(ds, MF_out, CO2_in, Ta_in, ps_in):
    """
    Purpose:
     Calculate CO2 mole fraction in uml/mol from mass density in mgCO2/m3.
    Usage:
     pfp_func_units.mgCO2pm3_to_umolpmol(ds, MF_out, CO2_in, Ta_in, ps_in)
    Author: PRI
    Date: August 2019
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [CO2_in, Ta_in, ps_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + MF_out + " not calculated"
            logger.error(msg)
            return 0
    CO2 = pfp_utils.GetVariable(ds, CO2_in)
    CO2 = pfp_utils.convert_units_func(ds, CO2, "mg/m^3")
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    Ta = pfp_utils.convert_units_func(ds, Ta, "degC")
    ps = pfp_utils.GetVariable(ds, ps_in)
    ps = pfp_utils.convert_units_func(ds, ps, "kPa")
    MF = pfp_utils.GetVariable(ds, MF_out)
    MF["Data"] = pfp_mf.co2_ppmfrommgCO2pm3(CO2["Data"], Ta["Data"], ps["Data"])
    MF["Flag"] = numpy.where(numpy.ma.getmaskarray(MF["Data"]) == True, ones, zeros)
    MF["Attr"]["units"] = "umol/mol"
    pfp_utils.CreateVariable(ds, MF)
    return 1

def mmolpm3_to_gH2Opm3(ds, AH_out, H2O_in):
    """
    Purpose:
     Function to convert mmol/m^3 (molar density) to g/m^3 (mass density).
    Usage:
     pfp_func_units.mmolpm3_to_gpm3(ds, AH_out, H2O_in)
    Author: PRI
    Date: August 2020
    """
    for item in [H2O_in]:
        if item not in list(ds.root["Variables"].keys()):
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
            return 0
    var_in = pfp_utils.GetVariable(ds, H2O_in)
    got_variance = False
    if var_in["Label"][-3:] == "_Vr" and var_in["Attr"]["units"] == "mmol^2/m^6":
        got_variance = True
        var_in["Data"] = numpy.ma.sqrt(var_in["Data"])
        var_in["Attr"]["units"] = "mmol/m^3"
    var_out = pfp_utils.convert_units_func(ds, var_in, "g/m^3", mode="quiet")
    var_out["Label"] = AH_out
    if got_variance:
        var_out["Data"] = var_out["Data"]*var_out["Data"]
        var_out["Attr"]["units"] = "g^2/m^6"
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def mmolpmol_to_gH2Opm3(ds, AH_out, MF_in, Ta_in, ps_in):
    """
    Purpose:
     Function to calculate absolute humidity given the water vapour mole
     fraction, air temperature and pressure.  Absolute humidity is not calculated
     if any of the input series are missing or if the specified output series
     already exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     pfp_func_units.mmolpmol_to_gpm3(ds,"AH_IRGA_Av","H2O_IRGA_Av","Ta_HMP_2m","ps")
    Author: PRI
    Date: September 2015
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [MF_in, Ta_in, ps_in]:
        if item not in list(ds.root["Variables"].keys()):
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
            return 0
    MF = pfp_utils.GetVariable(ds, MF_in)
    MF = pfp_utils.convert_units_func(ds, MF, "mmol/mol")
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    Ta = pfp_utils.convert_units_func(ds, Ta, "degC")
    ps = pfp_utils.GetVariable(ds, ps_in)
    ps = pfp_utils.convert_units_func(ds, ps, "kPa")
    AH = pfp_utils.GetVariable(ds, AH_out)
    AH["Data"] = pfp_mf.h2o_gpm3frommmolpmol(MF["Data"], Ta["Data"], ps["Data"])
    AH["Flag"] = numpy.where(numpy.ma.getmaskarray(AH["Data"]) == True, ones, zeros)
    AH["Attr"]["units"] = "g/m^3"
    pfp_utils.CreateVariable(ds, AH)
    return 1

def percent_to_mmolpmol(ds, MF_out, RH_in, Ta_in, ps_in):
    """
    Purpose:
     Calculate H2O mole fraction from relative humidity (RH).
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs,dtype=numpy.int32)
    ones = numpy.ones(nRecs,dtype=numpy.int32)
    for item in [RH_in, Ta_in, ps_in]:
        if item not in list(ds.root["Variables"].keys()):
            msg = " Requested series " + item + " not found, " + MF_out + " not calculated"
            logger.error(msg)
            return 0
    # get the relative humidity and check the units
    RH = pfp_utils.GetVariable(ds, RH_in)
    RH = pfp_utils.convert_units_func(ds, RH, "percent")
    # get the temperature and check the units
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    Ta = pfp_utils.convert_units_func(ds, Ta, "degC")
    # get the absoulte humidity
    AH_data = pfp_mf.absolutehumidityfromrelativehumidity(Ta["Data"], RH["Data"])
    # get the atmospheric pressure and check the units
    ps = pfp_utils.GetVariable(ds, ps_in)
    ps = pfp_utils.convert_units_func(ds, ps, "kPa")
    # get the output variable (created in pfp_ts.DoFunctions())
    MF = pfp_utils.GetVariable(ds, MF_out)
    # do the business
    MF["Data"] = pfp_mf.h2o_mmolpmolfromgpm3(AH_data, Ta["Data"], ps["Data"])
    MF["Flag"] = numpy.where(numpy.ma.getmaskarray(MF["Data"]) == True, ones, zeros)
    MF["Attr"]["units"] = "mmol/mol"
    # put the output variable back into the data structure
    pfp_utils.CreateVariable(ds, MF)
    return 1

def percent_to_gH2Opm3(ds, AH_out, RH_in, Ta_in):
    """
    Purpose:
     Function to calculate absolute humidity given relative humidity and
     air temperature.  Absolute humidity is not calculated if any of the
     input series are missing or if the specified output series already
     exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     pfp_func_units.percent_to_gpm3(ds,"AH_HMP_2m","RH_HMP_2m","Ta_HMP_2m")
    Author: PRI
    Date: September 2015
    """
    nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [RH_in, Ta_in]:
        if item not in ds.root["Variables"].keys():
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
            return 0
    # get the relative humidity and check the units
    RH = pfp_utils.GetVariable(ds, RH_in)
    RH = pfp_utils.convert_units_func(ds, RH, "percent")
    # get the temperature and check the units
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    Ta = pfp_utils.convert_units_func(ds, Ta, "degC")
    # get the absolute humidity
    AH = pfp_utils.GetVariable(ds, AH_out)
    AH["Data"] = pfp_mf.absolutehumidityfromrelativehumidity(Ta["Data"], RH["Data"])
    AH["Flag"] = numpy.where(numpy.ma.getmaskarray(AH["Data"]) == True, ones, zeros)
    AH["Attr"]["units"] = "g/m^3"
    pfp_utils.CreateVariable(ds, AH)
    return 1

def Pa_to_kPa(ds, ps_out, ps_in):
    """
    Purpose:
     Function to convert pressure from Pa to kPa.
    Usage:
     pfp_func_units.Pa_to_kPa(ds, ps_out, ps_in)
    Author: PRI
    Date: February 2018
    """
    var_in = pfp_utils.GetVariable(ds, ps_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "kPa", mode="quiet")
    var_out["Label"] = ps_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def percent_to_m3pm3(ds, Sws_out, Sws_in):
    """
    Purpose:
     Function to convert Sws in units of "percent" (1 to 100) to "frac" (0 to 1).
    Usage:
     pfp_func_units.percent_to_m3pm3(ds, Sws_out, Sws_in)
    Author: PRI
    Date: April 2020
    """
    var_in = pfp_utils.GetVariable(ds, Sws_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "m^3/m^3", mode="quiet")
    var_out["Label"] = Sws_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

#def Linear(ds, label_out, label_in, slope, offset):
    #"""
    #Purpose:
     #Function to apply a linear correction to a variable.
    #Usage:
     #pfp_func_units.Linear(ds, label_out, label_in, slope, offset)
    #Author: PRI
    #Date: August 2019
    #"""
    #nRecs = int(ds.root["Attributes"]["nc_nrecs"])
    #zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    #ones = numpy.ones(nRecs, dtype=numpy.int32)
    #for item in [label_in]:
        #if item not in ds.root["Variables"].keys():
            #msg = " Requested series " + item + " not found, " + label_out + " not calculated"
            #logger.error(msg)
            #return 0
    #var_in = pfp_utils.GetVariable(ds, label_in)
    #var_out = pfp_utils.GetVariable(ds, label_out)
    #var_out["Data"] = var_in["Data"] * float(slope) + float(offset)
    #var_out["Flag"] = numpy.where(numpy.ma.getmaskarray(var_out["Data"]) == True, ones, zeros)
    #pfp_utils.CreateVariable(ds, var_out)
    #return 1
