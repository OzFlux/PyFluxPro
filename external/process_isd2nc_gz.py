#!/usr/bin/env python3
# Python modules
import calendar
from collections import OrderedDict
import copy
import datetime
import gzip
import logging
import os
import sys
# 3rd party modules
from configobj import ConfigObj
import numpy
import pandas
import pytz
import scipy
import scipy.stats
import xlrd
import xlwt
# check the scripts directory is present
#if not os.path.exists("scripts"):
#    print("era52nc: the scripts directory is missing")
#    sys.exit()
# since the scripts directory is there, try importing the modules
#sys.path.append('scripts')
# pfp modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("external")-1])
import scripts.constants as c
import scripts.meteorologicalfunctions as mf
import scripts.pfp_ck as pfp_ck
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_ts as pfp_ts
import scripts.pfp_utils as pfp_utils
import scripts.pfp_compliance as pfp_compliance

now = datetime.datetime.now()
log_file_name = "isd2nc_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_path = "logfiles"
if not os.path.isdir(log_file_path):
    os.mkdir(log_file_path)
log_file_uri = os.path.join(log_file_path, log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_uri, to_screen=True)

def read_isd_file_gz(isd_file_path):
    """
    Purpose:
     Reads an ISD CSV file (gz or uncompressed) and returns the data in a data structure.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    msg = "Reading ISD file "+isd_file_path
    logger.info(msg)
    isd_file_name = os.path.split(isd_file_path)[1]
    isd_site_id = isd_file_name.split("-")
    isd_site_id = isd_site_id[0]+"-"+isd_site_id[1]
    with gzip.open(isd_file_path, 'rb') as fp:
            content = fp.readlines()
    # get a data structure
    ds = pfp_io.DataStructure()
    # get the site latitude, longitude and altitude
    ds.root["Attributes"]["altitude"] = float(content[0][46:51].decode('utf-8'))
    ds.root["Attributes"]["latitude"] = float(content[0][28:34].decode('utf-8'))/float(1000)
    ds.root["Attributes"]["longitude"] = float(content[0][34:41].decode('utf-8'))/float(1000)
    ds.root["Attributes"]["isd_site_id"] = isd_site_id
    ds.root["Attributes"]["processing_level"] = "L1"
    # initialise the data structure
    ds.root["Variables"]["DateTime"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Datetime","units":"none","statistic_type":"none"}}
    ds.root["Variables"]["Wd"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Wind direction","units":"degrees","statistic_type":"average"}}
    ds.root["Variables"]["Ws"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Wind speed","units":"m/s","statistic_type":"average"}}
    ds.root["Variables"]["Ta"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Air temperature","units":"degC","statistic_type":"average"}}
    ds.root["Variables"]["Td"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Dew point temperature","units":"degC","statistic_type":"average"}}
    ds.root["Variables"]["ps"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Surface pressure","units":"kPa","statistic_type":"average"}}
    ds.root["Variables"]["Precip"] = {"Data":[],"Flag":[],"Attr":{"long_name":"Precipitation","units":"mm","statistic_type":"sum"}}
    # define the codes for good data in the ISD file
    OK_obs_code = ["AUTO ","CRN05","CRN15","FM-12","FM-15","FM-16","SY-MT"]
    # iterate over the lines in the file and decode the data
    for i in range(len(content)-1):
    #for i in range(10):
        # filter out anything other than hourly data
        if content[i][41:46].decode('utf-8') not in OK_obs_code: continue
        YY = int(content[i][15:19].decode('utf-8'))
        MM = int(content[i][19:21].decode('utf-8'))
        DD = int(content[i][21:23].decode('utf-8'))
        HH = int(content[i][23:25].decode('utf-8'))
        mm = int(content[i][25:27].decode('utf-8'))
        dt = datetime.datetime(YY,MM,DD,HH,mm,0)
        ds.root["Variables"]["DateTime"]["Data"].append(pytz.utc.localize(dt))
        # wind direction, degT
        try:
            ds.root["Variables"]["Wd"]["Data"].append(float(content[i][60:63].decode('utf-8')))
        except:
            ds.root["Variables"]["Wd"]["Data"].append(float(999))
        # wind speed, m/s
        try:
            ds.root["Variables"]["Ws"]["Data"].append(float(content[i][65:69].decode('utf-8'))/float(10))
        except:
            ds.root["Variables"]["Ws"]["Data"].append(float(999.9))
        # air temperature, C
        try:
            ds.root["Variables"]["Ta"]["Data"].append(float(content[i][87:92].decode('utf-8'))/float(10))
        except:
            ds.root["Variables"]["Ta"]["Data"].append(float(999.9))
        # dew point temperature, C
        try:
            ds.root["Variables"]["Td"]["Data"].append(float(content[i][93:98].decode('utf-8'))/float(10))
        except:
            ds.root["Variables"]["Td"]["Data"].append(float(999.9))
        # sea level pressure, hPa
        try:
            ds.root["Variables"]["ps"]["Data"].append(float(content[i][99:104].decode('utf-8'))/float(10))
        except:
            ds.root["Variables"]["ps"]["Data"].append(float(9999.9))
        # precipitation, mm
        if content[i][108:111].decode('utf-8') == "AA1":
            try:
                ds.root["Variables"]["Precip"]["Data"].append(float(content[i][113:117].decode('utf-8'))/float(10))
            except:
                ds.root["Variables"]["Precip"]["Data"].append(float(999.9))
        else:
            ds.root["Variables"]["Precip"]["Data"].append(float(999.9))

    # add the time zone to the DateTime ataributes
    ds.root["Variables"]["DateTime"]["Attr"]["time_zone"] = "UTC"
    # convert from lists to masked arrays
    f0 = numpy.zeros(len(ds.root["Variables"]["DateTime"]["Data"]))
    f1 = numpy.ones(len(ds.root["Variables"]["DateTime"]["Data"]))
    ds.root["Variables"]["DateTime"]["Data"] = numpy.array(ds.root["Variables"]["DateTime"]["Data"])
    ds.root["Variables"]["DateTime"]["Flag"] = f0
    ds.root["Attributes"]["nc_nrecs"] = len(f0)
    dt_delta = pfp_utils.get_timestep(ds)
    ts = scipy.stats.mode(dt_delta)[0]/60
    ds.root["Attributes"]["time_step"] = ts[0]

    ds.root["Variables"]["Wd"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["Wd"]["Data"],999)
    ds.root["Variables"]["Wd"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["Wd"]["Data"])==True,f1,f0)
    ds.root["Variables"]["Ws"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["Ws"]["Data"],999.9)
    ds.root["Variables"]["Ws"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["Ws"]["Data"])==True,f1,f0)
    ds.root["Variables"]["Ta"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["Ta"]["Data"],999.9)
    ds.root["Variables"]["Ta"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["Ta"]["Data"])==True,f1,f0)
    ds.root["Variables"]["Td"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["Td"]["Data"],999.9)
    ds.root["Variables"]["Td"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["Td"]["Data"])==True,f1,f0)
    # hPa to kPa
    ds.root["Variables"]["ps"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["ps"]["Data"],9999.9)/float(10)
    ds.root["Variables"]["ps"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["ps"]["Data"])==True,f1,f0)
    # convert sea level pressure to station pressure
    site_altitude = float(ds.root["Attributes"]["altitude"])
    cfac = numpy.ma.exp((-1*site_altitude)/((ds.root["Variables"]["Ta"]["Data"]+273.15)*29.263))
    ds.root["Variables"]["ps"]["Data"] = ds.root["Variables"]["ps"]["Data"]*cfac
    # do precipitation and apply crude limits
    ds.root["Variables"]["Precip"]["Data"] = numpy.ma.masked_equal(ds.root["Variables"]["Precip"]["Data"],999.9)
    condition = (ds.root["Variables"]["Precip"]["Data"]<0)|(ds.root["Variables"]["Precip"]["Data"]>100)
    ds.root["Variables"]["Precip"]["Data"] = numpy.ma.masked_where(condition,ds.root["Variables"]["Precip"]["Data"])
    ds.root["Variables"]["Precip"]["Flag"] = numpy.where(numpy.ma.getmaskarray(ds.root["Variables"]["Precip"]["Data"])==True,f1,f0)
    # get the humidities from Td
    pfp_ts.RelativeHumidityFromDewpoint(ds)
    pfp_ts.AbsoluteHumidityFromRelativeHumidity(ds)
    pfp_ts.SpecificHumidityFromAbsoluteHumidity(ds)

    # return the data
    return ds

def read_isd_file_csv(isd_file_path):
    """
    Purpose:
     Reads a NOAA ISD CSV file downlaoded from https://www.ncei.noaa.gov/data/global-hourly/access/
     These files used to be field formatted ASCII where the character position in a line
     of ASCII determined the data type.  Some time in 2020 or 2021, the old FFA format
     was replaced with CSV.
     The format of the old-style .gz files is described in
     https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
     This document still describes the data in the new CSV format.
    Usage:
    Side effects:
     Returns a PFP data structure with the data at the site time step.
    Author: PRI
    Date: July 2021
    """
    msg = " Reading " + isd_file_path
    logger.info(msg)
    # list of variables to read from the CSV file
    csv_labels = ["STATION", "DATE", "LATITUDE", "LONGITUDE", "ELEVATION", "REPORT_TYPE",
                   "QUALITY_CONTROL", "WND", "TMP", "DEW", "SLP", "AA1", "AA2", "AA3", "AA4"]
    # read the CSV file
    df = pandas.read_csv(isd_file_path, delimiter=",", header=0)
    # remove items from csv_labels that are not in the data frame
    df_labels = df.columns.to_list()
    for csv_label in list(csv_labels):
        if csv_label not in df_labels:
            csv_labels.remove(csv_label)
    # keep only what we need
    df = df[csv_labels]
    # remove duplicate dates, keep the SYNOP (FM-12) reports
    # first, we find the duplicate dates
    df["Duplicates"] = df["DATE"].duplicated()
    # next, we drop rows with duplicate dates that are not SYNOP reports
    df = df.drop(df[(df["Duplicates"]) & (df["REPORT_TYPE"] != "FM-12")].index)
    # then check for duplicates again
    df["Duplicates"] = df["DATE"].duplicated()
    if df["Duplicates"].sum() != 0:
        msg = " Unable to remove all duplicate dates in files"
        logger.error(msg)
        raise ValueError(msg)
    # convert the date in the CSV file to a pandas datetime
    df["TIMESTAMP"] = pandas.to_datetime(df["DATE"].astype("string"),errors="raise")
    # find all of the timestamps (should only be 1)
    timestamps = list(df.select_dtypes(include=['datetime64']))
    # take the first if more than 1
    timestamp = timestamps[0]
    # use the timestamp as the index
    df.set_index(timestamp, inplace=True)
    df.index = df.index.round('1S')
    # wind direction field, see isd_format_document.pdf for details
    wind = df["WND"].str.split(',', expand=True)
    df["Wd"] = wind[0].apply(pandas.to_numeric, errors='coerce')
    df["Ws"] = wind[3].apply(pandas.to_numeric, errors='coerce')/float(10)
    del df["WND"]
    # air temperature
    temperature = df["TMP"].str.split(',', expand=True)
    df["Ta"] = temperature[0].apply(pandas.to_numeric, errors='coerce')/float(10)
    del df["TMP"]
    # dew point temperature
    dew_point = df["DEW"].str.split(',', expand=True)
    df["Td"] = dew_point[0].apply(pandas.to_numeric, errors='coerce')/float(10)
    del df["DEW"]
    # surface pressure
    surface_pressure = df["SLP"].str.split(',', expand=True)
    df["ps"] = surface_pressure[0].apply(pandas.to_numeric, errors='coerce')/float(100)
    del df["SLP"]
    # Precipitation is stored in columns AA1 to AA4 but not all columns will be present
    #
    # Within each column, precipitation is stored as "HH,PPPP,C,Q" where HH is the
    # period over which the precipitation was accumulated (e.g. 1, 3, 6, 24 hours),
    # PPPP is the precipitation amount in mm*10, C is the condition code and Q is
    # the QC flag (1 = passed all QC checks)
    #
    # Column AA1 contains most of the precipitation data.  When precipitation data is
    # available for 2 accumulation periods e.g. 3 hours and 6 or 24 hours, the second
    # accumulation period is given in AA2.  And so on for up to 4 separate accumulation'
    # periods e.g. 1 hour, 3 hours, 6 hours and 24 hours.
    #
    # get a list of the precipitation columns in the data frame
    precip_labels = [l for l in df.columns.to_list() if "AA" in l]
    # create a data frame for the precipitation data, same index as main data frame
    df_precip = pandas.DataFrame(index=df.index)
    # loop over the precipitation fields
    for precip_label in precip_labels:
        # split the "HH,PPPP,C,Q" fields to get individual parts
        tmp = df[precip_label].str.split(',', expand=True)
        # name the columns
        tmp.columns = ["Period", "Amount", "Condition", "Quality"]
        # coerce to numeric values
        tmp = tmp.apply(pandas.to_numeric, errors='coerce')
        # loop over the accumulation periods
        for n in [1, 3, 6, 24]:
            # get the data for this accumulation period and store in a new column
            # e.g. "3_hourly_AA1"
            tmp.loc[(tmp["Period"] == n) & (tmp["Quality"] == 1),
                     str(n)+"_hourly_"+precip_label] = tmp["Amount"]
        # drop the intermediate columns, no longer needed
        tmp = tmp.drop(["Period", "Amount", "Condition", "Quality"], axis=1)
        # concatenate the new data
        df_precip = pandas.concat([df_precip, tmp], axis=1)
        # drop the individual columns e.g. AA1, AA2 etc
        df.drop(precip_label, axis=1, inplace=True)
    # now loop over the accumulation periods and combine to get a single column
    # for each accumulation period
    for n in [1, 3, 6, 24]:
        # list of column headings for this accumulation period
        label = str(n)+"_hourly"
        hour_labels = [l for l in df_precip.columns.to_list() if label in l]
        # rename the first column e.g. "3_hourly_AA1" to "3_hourly"
        df_precip.rename({hour_labels[0]: label}, axis=1, inplace=True)
        # loop over the remaining columns and merge into a single column for this
        # accumulation period
        for hour_label in hour_labels[1:]:
            # merge "3_hourly" with "3_hourly_AA2" etc
            df_precip[label] = df_precip[label].combine_first(df_precip[hour_label])
            # convert mm*10 to mm
            df_precip[label] = df_precip[label]/float(10)
            # delete columns that are no longer needed
            df_precip.drop(hour_label, axis=1, inplace=True)
    # print the sum of the 1, 3, 6 and 24 hourly accumulation periods (we expect them to
    # be equal)
    msg = " 1 hourly precipitation total is " + str(round(df_precip["1_hourly"].sum(), 4))
    logger.info(msg)
    msg = " 3 hourly precipitation total is " + str(round(df_precip["3_hourly"].sum(), 4))
    logger.info(msg)
    msg = " 6 hourly precipitation total is " + str(round(df_precip["6_hourly"].sum(), 4))
    logger.info(msg)
    msg = " 24 hourly precipitation total is " + str(round(df_precip["24_hourly"].sum(), 4))
    logger.info(msg)
    # choose the most common accumulation period
    msg = " Using " + df_precip.count().idxmax() + " for precipitation"
    logger.info(msg)
    # and use it for the precipitation data
    df["Precip"] = df_precip[df_precip.count().idxmax()]
    # now copy the data from a pandas data frame to a PFP data structure
    nrecs = len(df)
    ones = numpy.ones(nrecs)
    zeros = numpy.zeros(nrecs)
    # create a data structure
    ds_its = pfp_io.DataStructure()
    # set the global attributes
    ds_its.root["Attributes"]["nc_nrecs"] = nrecs
    ds_its.root["Attributes"]["altitude"] = float(df["ELEVATION"][0])
    ds_its.root["Attributes"]["latitude"] = float(df["LATITUDE"][0])
    ds_its.root["Attributes"]["longitude"] = float(df["LONGITUDE"][0])
    ds_its.root["Attributes"]["isd_site_id"] = int(df["STATION"][0])
    # get the datetime variable
    ldt = pfp_utils.CreateEmptyVariable("DateTime", nrecs)
    ldt["Data"] = numpy.array(df.index.to_pydatetime())
    ldt["Flag"] = zeros
    ldt["Attr"] = {"long_name": "Datetime in UTC", "units": ""}
    pfp_utils.CreateVariable(ds_its, ldt)
    # get the time step
    dt = pfp_utils.get_timestep(ds_its)
    time_step = int(scipy.stats.mode(dt/float(60))[0][0])
    if time_step not in [10, 30, 60, 180]:
        msg = " Time step (" + str(time_step) + ") must be 10, 30, 60 or 180 minutes"
        logger.error(msg)
        raise ValueError(msg)
    else:
        ds_its.root["Attributes"]["time_step"] = int(scipy.stats.mode(dt/float(60))[0][0])
    # now add the other variables
    # wind direction
    Wd = pfp_utils.CreateEmptyVariable("Wd", nrecs, datetime=ldt["Data"])
    Wd["Data"] = numpy.ma.masked_equal(df["Wd"].values, 999)
    Wd["Flag"] = numpy.where(numpy.ma.getmaskarray(Wd["Data"]) == True, ones, zeros)
    Wd["Attr"] = {"long_name": "Wind direction", "statistic_type": "average",
                  "standard_name": "wind_from_direction", "units": "degrees"}
    pfp_utils.CreateVariable(ds_its, Wd)
    # wind speed
    Ws = pfp_utils.CreateEmptyVariable("Ws", nrecs, datetime=ldt["Data"])
    Ws["Data"] = numpy.ma.masked_equal(df["Ws"].values, 999.9)
    Ws["Flag"] = numpy.where(numpy.ma.getmaskarray(Ws["Data"]) == True, ones, zeros)
    Ws["Attr"] = {"long_name": "Wind speed", "statistic_type": "average",
                  "standard_name": "wind_speed", "units": "m/s"}
    pfp_utils.CreateVariable(ds_its, Ws)
    # air temperature
    Ta = pfp_utils.CreateEmptyVariable("Ta", nrecs, datetime=ldt["Data"])
    Ta["Data"] = numpy.ma.masked_equal(df["Ta"].values, 999.9)
    Ta["Flag"] = numpy.where(numpy.ma.getmaskarray(Ta["Data"]) == True, ones, zeros)
    Ta["Attr"] = {"long_name": "Air temperature", "statistic_type": "average",
                  "standard_name": "air_temperature", "units": "degC"}
    pfp_utils.CreateVariable(ds_its, Ta)
    # dew point temperature
    Td = pfp_utils.CreateEmptyVariable("Td", nrecs, datetime=ldt["Data"])
    Td["Data"] = numpy.ma.masked_equal(df["Td"].values, 999.9)
    Td["Flag"] = numpy.where(numpy.ma.getmaskarray(Td["Data"]) == True, ones, zeros)
    Td["Attr"] = {"long_name": "Dew point temperature", "statistic_type": "average",
                  "standard_name": "dew_point_temperature", "units": "degC"}
    pfp_utils.CreateVariable(ds_its, Td)
    # surface pressure
    ps = pfp_utils.CreateEmptyVariable("ps", nrecs, datetime=ldt["Data"])
    site_altitude = float(ds_its.root["Attributes"]["altitude"])
    cfac = numpy.ma.exp((-1*site_altitude)/((Ta["Data"]+273.15)*29.263))
    ps["Data"] = numpy.ma.masked_equal(df["ps"].values, 9999.9)
    ps["Data"] = ps["Data"]*cfac
    ps["Flag"] = numpy.where(numpy.ma.getmaskarray(ps["Data"]) == True, ones, zeros)
    ps["Attr"] = {"long_name": "Surface pressure", "statistic_type": "average",
                  "standard_name": "surface_air_pressure", "units": "kPa"}
    pfp_utils.CreateVariable(ds_its, ps)
    # precipitation
    Precip = pfp_utils.CreateEmptyVariable("Precip", nrecs, datetime=ldt["Data"])
    Precip["Data"] = numpy.ma.masked_equal(df["Precip"].values, 999.9)
    Precip["Flag"] = numpy.where(numpy.ma.getmaskarray(Precip["Data"]) == True, ones, zeros)
    Precip["Attr"] = {"long_name": "Rainfall", "statistic_type": "sum",
                      "standard_name": "thickness_of_rainfall_amount", "units": "mm"}
    pfp_utils.CreateVariable(ds_its, Precip)
    # get the humidities from Td
    pfp_ts.RelativeHumidityFromDewpoint(ds)
    pfp_ts.AbsoluteHumidityFromRelativeHumidity(ds)
    pfp_ts.SpecificHumidityFromAbsoluteHumidity(ds)

    # return the data
    return ds_its

def interpolate_ds(ds_in, ts, k=3):
    """
    Purpose:
     Interpolate the contents of a data structure onto a different time step.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    # instance the output data structure
    ds_out = pfp_io.DataStructure()
    # copy the global attributes
    for key in list(ds_in.root["Attributes"].keys()):
        ds_out.root["Attributes"][key] = ds_in.root["Attributes"][key]
    # add the time step
    ds_out.root["Attributes"]["time_step"] = str(ts)
    # generate a regular time series at the required time step
    dt = ds_in.root["Variables"]["DateTime"]["Data"]
    dt0 = dt[0] - datetime.timedelta(minutes=30)
    start = datetime.datetime(dt0.year, dt0.month, dt0.day, dt0.hour, 0, 0)
    dt1 = dt[-1] + datetime.timedelta(minutes=30)
    end = datetime.datetime(dt1.year, dt1.month, dt1.day, dt1.hour, 0, 0)
    idt = [result for result in perdelta(start, end, datetime.timedelta(minutes=ts))]
    x1 = numpy.array([toTimestamp(dt[i]) for i in range(len(dt))])
    x2 = numpy.array([toTimestamp(idt[i]) for i in range(len(idt))])
    # loop over the series in the data structure and interpolate
    ds_out.root["Variables"]["DateTime"] = {}
    ds_out.root["Variables"]["DateTime"]["Data"] = idt
    ds_out.root["Variables"]["DateTime"]["Flag"] = numpy.zeros(len(idt))
    ds_out.root["Variables"]["DateTime"]["Attr"] = {"long_name":"Datetime","units":"none"}
    ds_out.root["Attributes"]["nc_nrecs"] = len(idt)
    series_list = list(ds_in.root["Variables"].keys())
    if "DateTime" in series_list:
        series_list.remove("DateTime")
    for label in series_list:
        #print label
        data_in = pfp_utils.GetVariable(ds_in, label)
        # check if we are dealing with precipitation
        if "Precip" in label:
            # precipitation shouldn't be interpolated, just assign any precipitation
            # to the ISD time stamp.
            data_out = numpy.ma.zeros(len(idt), dtype=numpy.float64)
            idx = numpy.searchsorted(x2, numpy.intersect1d(x2, x1))
            idy = numpy.searchsorted(x1, numpy.intersect1d(x1, x2))
            data_out[idx] = data_in["Data"][idy]
        else:
            # interpolate everything else
            data_out = interpolate_1d(x1, data_in["Data"], x2)
        flag_out = numpy.zeros(len(idt))
        attr_out = attr_in
        # pfp_utils.CreateSeries(ds_out, label, data_out, Flag=flag_out, Attr=attr_out)
        pfp_utils.CreateVariable(ds_out, label)
    return ds_out

def interpolate_1d(x1, y1, x2):
    """
    Purpose:
     Interpolate data from one time step to another.
    Assumptions:
    Usage:
    Author: PRI
    Date: June 2017
    """
    # off we go
    if numpy.ma.is_masked(y1):
        # check we have at least 2 non-masked points
        if numpy.ma.count(y1) >= 2:
            # input Y array is a masked array
            idx = numpy.where(numpy.ma.getmaskarray(y1) == False)[0]
            int_fn = scipy.interpolate.Akima1DInterpolator(x1[idx], y1[idx].data)
            y2 = int_fn(x2)
        else:
            msg = "Not enough points (<2) to interpolate"
            logger.warning(msg)
            y2 = numpy.ma.ones(len(x2))*float(c.missing_value)
    else:
        int_fn = scipy.interpolate.Akima1DInterpolator(x1, y1)
        y2 = int_fn(x2)
    return y2

def perdelta(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def toTimestamp(d):
    return calendar.timegm(d.timetuple())

def convert_time_zone(ds, from_time_zone, to_time_zone):
    """
    Purpose:
     Convert the datetime series in a data structure from one timezone to another.
    Usage:
    Author: PRI
    Date: June 2017
    """
    local_time_zone = pytz.timezone(to_time_zone)
    from_dt_naive = ds.root["Variables"]["DateTime"]["Data"]
    to_dt_aware = [dt.replace(tzinfo=from_time_zone).astimezone(local_time_zone) for dt in from_dt_naive]
    to_dt_naive = [dt.replace(tzinfo=None) for dt in to_dt_aware]
    ds.root["Variables"]["DateTime"]["Data"] = to_dt_naive
    ds.root["Variables"]["DateTime"]["Attr"]["time_zone"] = local_time_zone
    return

def read_site_master(xl_file_path, sheet_name):
    """
    """
    xl_book = xlrd.open_workbook(xl_file_path)
    xl_sheet = xl_book.sheet_by_name(sheet_name)
    last_row = int(xl_sheet.nrows)
    # find the header and first data rows
    for i in range(last_row):
        if xl_sheet.cell(i,0).value == "Site":
            header_row = i
            first_data_row = header_row + 1
            break
    # read the header row
    header_row_values = xl_sheet.row_values(header_row)
    # read the site data from the master Excel spreadsheet
    site_info = OrderedDict()
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = OrderedDict()
        site_info[site_name]["site_name"] = site_name
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value

    return site_info

def remove_duplicates(ds):
    """
    Remove duplicate timestamps, similar to Peter's solution in pfp_ts.py
    MergeDataStructures at L1
    """
    # get the datetime
    dtn = pfp_utils.GetVariable(ds, "DateTime")
    # remove duplicate timestamps
    dtn_unique, index_unique = numpy.unique(dtn["Data"], return_index=True)
    # restore the original order of the unique timestamps
    dtn_sorted = dtn_unique[numpy.argsort(index_unique)]
    # check to see if there were duplicates
    if len(dtn_sorted) < len(dtn["Data"]):
        n = len(dtn["Data"]) - len(dtn_sorted)
        msg = str(n) + " duplicate time stamps were removed for isd site "
        logger.warning(msg)
    nrecs = len(dtn_sorted)
    labels = list(ds.root["Variables"].keys())
    #if "DateTime" in labels:
    #    labels.remove("DateTime")
    for label in labels:
        var1 = pfp_utils.CreateEmptyVariable(label, nrecs)
        varn = pfp_utils.GetVariable(ds, label)
        var1["Data"] = varn["Data"][index_unique]
        var1["Flag"] = varn["Flag"][index_unique]
        var1["Attr"] = varn["Attr"]
        pfp_utils.CreateVariable(ds, var1)
    return ds

def xl_write_ISD_timesteps(xl_file_path, data):
    """
    Purpose:
     Writes a dictionary to a worksheet in an Excel workbook.
     This routine has 2 arguments,an Excel worksheet instance and
     a dictionary of data to be written out.  The dictionary
     format needs to be:
      data[site][year]["mean"]
      data[site][year]["stdev"]
      data[site][year]["mode"]
    Usage:
     pfp_io.xl_write_ISD_timesteps(xl_file_path, data)
      where xl_file_path is an Excel workbook file name
            data         is a dictionary as defined above
    Side effects:
     Writes to an Excel worksheet instance.
    Called by:
    Calls:
    Author: PRI
    Date: August 2017
    """
    # get a workbook
    xl_book = xlwt.Workbook()
    # get a list of the sheets to add
    site_list = list(data.keys())
    year_list = list(data[site_list[0]].keys())
    stat_list = list(data[site_list[0]][year_list[0]].keys())
    # loop over the statistics
    for stat in stat_list:
        # add a worksheet for the statistics
        xl_sheet = xl_book.add_sheet(stat)
        # write the header line
        for col, year in enumerate(year_list):
            xl_sheet.write(0,col+1,year)
        # write the data, one row per site, one column per year
        for row, site in enumerate(site_list):
            xl_sheet.write(row+1,0,site)
            for col, year in enumerate(year_list):
                if stat in list(data[site][year].keys()):
                    xl_sheet.write(row+1,col+1,data[site][year][stat])
                else:
                    xl_sheet.write(row+1,col+1,"")
    # save the workbook
    xl_book.save(xl_file_path)

    return

# =============================================== MAIN =========================================

# read the control file file
cfg_file_path = "process_isd2nc.txt"
msg = " Loading the control file"
logger.info(msg)
cfg = pfp_io.get_controlfilecontents(cfg_file_path, mode="verbose")
xl_file_path = cfg["Files"]["xl_file_path"]
xl_sheet_name = cfg["Files"]["xl_sheet_name"]
isd_base_path = cfg["Files"]["isd_base_path"]
out_base_path = cfg["Files"]["out_base_path"]
# read the site master spreadsheet
site_info = read_site_master(xl_file_path, xl_sheet_name)
# get a list of sites
site_list = list(site_info.keys())
# creat a dictionary to hold the ISD site time steps
isd_time_steps = OrderedDict()

for site in site_list:
    # construct the output file path
        nc_out_path = os.path.join(out_base_path,site,"Data","ISD",site+"_ISD.nc")
    # construct the config dictionary for the concatenate routine
    cfg_concat = ConfigObj(indent_type="    ")
    cfg_concat["Options"] = {"NumberOfDimensions":1,
                            "MaxGapInterpolate":0,
                            "FixTimeStepMethod":"round",
                            "Truncate":"No",
                            "TruncateThreshold":50,
                            "SeriesToCheck":[]}
    cfg_concat["Files"] = {"Out":{"ncFileName":nc_out_path},
                          "In":{}}
    # get the list of ISD stations to be used for this site
    isd_site_list = []
    for item in ["ISD_ID_1","ISD_ID_2","ISD_ID_3","ISD_ID_4"]:
        if len(site_info[site][item]) != 0:
            isd_site_list.append(site_info[site][item])

    # now get a dictionary that ties the ISD station ID to a number
    # that will be appended to the variable name
    site_index = {}
    for n, isd_site in enumerate(isd_site_list):
        site_index[isd_site] = n
    if not isinstance(isd_site_list, list):
        isd_site_list = [isd_site_list]
    time_zone = site_info[site]["Time zone"]
    time_step = int(round(float(site_info[site]["Time step"])))
    start_year = int(site_info[site]["Start year"])
    try:
        end_year = int(site_info[site]["End year"])
    except:
        end_year = int(now.year)
    # get the list of years to process
    year_list = list(range(start_year,end_year+1))
    for n, year in enumerate(year_list):
        # we will collect the data for each site for this year into a single dictionary
        ds_out = {}
        isd_year_path = os.path.join(isd_base_path,str(year))
        # check if files in isd_year_path
        file_year_list = os.listdir(isd_year_path)
        print(file_year_list)
        if len(file_year_list) == 0:
            msg = " No files found in " + isd_year_path
            logger.error(msg)
            continue
        for isd_site in isd_site_list:
            if isd_site not in list(isd_time_steps.keys()):
                isd_time_steps[isd_site] = OrderedDict()
            if year not in list(isd_time_steps[isd_site].keys()):
                isd_time_steps[isd_site][year] = OrderedDict()
            # check if downloaded files are gz or csv files
            isd_site_year = str(isd_site)+"-"+str(year)
            if isd_site_year+".gz" in file_year_list:
            isd_file_path = os.path.join(isd_year_path,str(isd_site)+"-"+str(year)+".gz")
                ds_isd = read_isd_file_gz(isd_file_path)
            elif isd_site_year+".csv" in file_year_list:
                isd_file_path_csv = os.path.join(isd_year_path, isd_site.replace("-", "") + ".csv")
                ds_isd = read_isd_file_csv(isd_file_path)
            else:
                msg = " No gz or csv files found in " + isd_year_path
                logger.error(msg)
            ds_in = remove_duplicates(ds_isd)
            # get an array of time steps in seconds
            dt = pfp_utils.get_timestep(ds_in)
            # and get dt in minutes
            dt = dt/float(60)
            isd_time_steps[isd_site][year]["mean"] = numpy.mean(dt)
            isd_time_steps[isd_site][year]["stdev"] = numpy.std(dt)
            isd_time_steps[isd_site][year]["mode"] = scipy.stats.mode(dt)[0][0]

            # interpolate from the ISD site time step to the tower time step
            labels = [l for l in list(ds_in.root["Variables"].keys()) if "DateTime" not in l]
            # interpolate from the ACCESS time step (60 minutes) to the tower time step
            ds_out[site_index[isd_site]] = pfp_ts.InterpolateDataStructure(ds_in, labels=labels,
                                                            new_time_step=time_step,
                                                            sums="interpolate",
                                                            mode="quiet")
            
            # adjust time from UTC to local using the time zone
            convert_time_zone(ds_out[site_index[isd_site]], pytz.utc, time_zone)
            # add some useful global attributes
            ds_out[site_index[isd_site]].root["Attributes"]["isd_site_id"] = isd_site
            ds_out[site_index[isd_site]].root["Attributes"]["time_zone"] = time_zone
            # write out a netCDF file for each ISD site and each year
            #nc_file_name = isd_site+"_"+str(year)+".nc"
            #nc_dir_path = os.path.join(out_base_path,site,"Data","ISD")
            #if not os.path.exists(nc_dir_path):
                #os.makedirs(nc_dir_path)
            #nc_file_path = os.path.join(nc_dir_path,nc_file_name)
            #nc_file = pfp_io.nc_open_write(nc_file_path)
            #pfp_io.nc_write_series(nc_file, ds_out[site_index[isd_site]], ndims=1)
        # now we merge the data structures for each ISD station into a single data structure
        # first, instance a data structure
        ds_all = pfp_io.DataStructure()
        ds_all.root["Attributes"]["latitude"] = site_info[site]["Latitude"]
        ds_all.root["Attributes"]["longitude"] = site_info[site]["Longitude"]
        ds_all.root["Attributes"]["altitude"] = site_info[site]["Altitude"]
        ds_all.root["Attributes"]["site_name"] = site_info[site]["site_name"]
        # now loop over the data structures for each ISD station and get the earliest
        # start time and the latest end time
        start_datetime = []
        end_datetime = []
        for i in list(ds_out.keys()):
            # print(i)
            start_datetime.append(ds_out[i].root["Variables"]["DateTime"]["Data"][0])
            end_datetime.append(ds_out[i].root["Variables"]["DateTime"]["Data"][-1])
        start = min(start_datetime)
        end = max(end_datetime)
        # now make a datetime series at the required time step from the earliest start
        # datetime to the latest end datetime
        ldt_all = [result for result in perdelta(start, end, datetime.timedelta(minutes=time_step))]
        nrecs = len(ldt_all)
        ds_all.root["Attributes"]["nc_nrecs"] = nrecs
        # and add the datetime to the all-stations data structure
        ds_all.root["Variables"]["DateTime"] = {}
        ds_all.root["Variables"]["DateTime"]["Data"] = ldt_all
        ds_all.root["Variables"]["DateTime"]["Flag"] = numpy.zeros(len(ldt_all))
        ds_all.root["Variables"]["DateTime"]["Attr"] = {"long_name":"DateTime","units":"none"}
        # now copy the contents of the ISD station data structures to the all-stations
        # data structure
        for i in list(ds_out.keys()):
            isd_site_id = ds_out[i].root["Attributes"]["isd_site_id"]
            # copy the global attributes
            gattr_list = list(ds_out[i].root["Attributes"].keys())
            # remove the site specific global attributes
            for item in ["latitude", "longitude", "altitude", "isd_site_id"]:
                gattr_list.remove(item)
            # copy everything else
            for gattr in gattr_list:
                if gattr not in ds_all.root["Attributes"]:
                    ds_all.root["Attributes"][gattr] = ds_out[i].root["Attributes"][gattr]
            # and update the site specific global attributes
            ds_all.root["Attributes"]["time_zone_"+isd_site_id] = ds_out[i].root["Attributes"]["time_zone"]
            ds_all.root["Attributes"]["latitude_"+isd_site_id] = ds_out[i].root["Attributes"]["latitude"]
            ds_all.root["Attributes"]["longitude_"+isd_site_id] = ds_out[i].root["Attributes"]["longitude"]
            ds_all.root["Attributes"]["altitude_"+isd_site_id] = ds_out[i].root["Attributes"]["altitude"]
            # now copy the variables
            # first, we get the indices of matching datetimes
            ldt_one = ds_out[i].root["Variables"]["DateTime"]["Data"]
            idx = numpy.searchsorted(ldt_all, numpy.intersect1d(ldt_all, ldt_one))
            idy = numpy.searchsorted(ldt_one, numpy.intersect1d(ldt_one, ldt_all))
            # then we get a list of the variables to copy
            series_list = list(ds_out[i].root["Variables"].keys())
            # and remove the datetime
            if "DateTime" in series_list:
                series_list.remove("DateTime")
            # and then we loop over the variables to be copied
            for label in series_list:
                # append a number, unique to each ISD station, to the variable label
                all_label = label+"_"+str(i)
                # create empty data and flag arrays
                variable = pfp_utils.CreateEmptyVariable(all_label, nrecs)
                pfp_utils.CreateVariable(ds_all, variable)
                # read the data out of the ISD site data structure
                data = pfp_utils.GetVariable(ds_out[i], label)
                # add the ISD site ID
                attr = data["Attr"]
                attr["isd_site_id"] = isd_site_id
                # put the data, flag and attributes into the all-in-one data structure
                ds_all.root["Variables"][all_label]["Data"][idx] = data["Data"][idy]
                ds_all.root["Variables"][all_label]["Flag"][idx] = data["Flag"][idy]
                ds_all.root["Variables"][all_label]["Attr"] = copy.deepcopy(attr)
        # do the QC checks
        cfg_qc = copy.deepcopy(cfg)
        cfg_labels = list(cfg["Variables"].keys())
        for cfg_label in cfg_labels:
            ds_labels = [l for l in list(ds_all.root["Variables"].keys()) if l.split("_")[0] == cfg_label]
            for ds_label in ds_labels:
                cfg_qc["Variables"][ds_label] = copy.deepcopy(cfg["Variables"][cfg_label])
            cfg_qc["Variables"].pop(cfg_label)
        pfp_ck.do_qcchecks(cfg_qc, ds_all)
        # quick and dirty hack needed because V3.2 and V3.3 treat valid_range
        # attribute differently.
        labels = list(ds_all.root["Variables"].keys())
        for label in labels:
            attrs = list(ds_all.root["Variables"][label]["Attr"].keys())
            for attr in ["valid_range"]:
                if attr in attrs:
                    ds_all.root["Variables"][label]["Attr"].pop(attr)
        # write the netCDF file with the combined data for this year
            nc_dir_path = os.path.join(out_base_path,site,"Data","ISD")
            nc_file_name = site+"_ISD_"+str(year)+".nc"
        nc_file_path = os.path.join(nc_dir_path,nc_file_name)
        pfp_io.NetCDFWrite(nc_file_path, ds_all)
        cfg_concat["Files"]["In"][str(n)] = nc_file_path
    # concatenate the yearly files for this site
    info = pfp_compliance.ParseConcatenateControlFile(cfg_concat)
    pfp_io.NetCDFConcatenate(info)

# write the time steps out to an Excel file
xl_file_path = os.path.join(isd_base_path, "ISD_site_timesteps.xls")
xl_write_ISD_timesteps(xl_file_path, isd_time_steps)

logger.info("All done")
