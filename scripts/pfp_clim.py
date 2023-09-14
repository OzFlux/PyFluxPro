# standard modules
import calendar
import datetime
import logging
import os
# 3rd party modules
from scipy.interpolate import griddata
import numpy
import xlwt
# PFP modules
from scripts import constants as c
from scripts import pfp_io
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def do_2dinterpolation(array_2d, tile="no", method="linear"):
    """
    Takes a 2d array as input and;
     1) tiles this into a 3 x 3 space (9 repeats of the original 2d array in 3 columns and 3 rows)
     2) removes the missing data (c.missing_value) from the tiled array
     3) does a bi-linear interpolation to replace the the missing data
     4) returns the central tile
     The effect is to replace missing data in the original 2d array with data from a bi-linear
     interpolation, the tiling repeats the original array along its boundaries to avoid problems
     at the array edges.
     Ian McHugh's cleaned up version (avoids error messages from matplotlib version of griddata).
     Checked by comparing the Fci(day) values from the original code and from this version.  The
     values were the same so this version pushed to GitHub on 11/5/2015.
    """
    WasMA = False
    if numpy.ma.isMA(array_2d):
        WasMA = True
        array_2d = numpy.ma.filled(array_2d, float(c.missing_value))
    if tile.lower() == "yes":
        # Tile the 2d array into a 3 by 3 array
        array_2d = numpy.tile(array_2d, (3, 3))
    # Get the dimensions of the tiled array and create coordinates and grid
    num_x = numpy.shape(array_2d)[1]
    array_x = numpy.arange(0, num_x)
    num_y = numpy.shape(array_2d)[0]
    array_y = numpy.arange(0, num_y)
    coords_x, coords_y = numpy.meshgrid(array_x, array_y)
    # Make a flat array of the tiled data
    data_1d = array_2d.flatten()
    # Make a 2d array of the coordinates
    data_coords = numpy.column_stack([coords_x.flatten(), coords_y.flatten()])
    # Define an index that will return all valid data for the array
    index = numpy.where(data_1d != c.missing_value)
    # Do the interpolation
    try:
        grid_z = griddata(data_coords[index], data_1d[index], (coords_x, coords_y), method = 'linear')
        if tile.lower() == "yes":
            # Retrieve the central tile
            array_2d_filled = grid_z[num_y // 3: num_y // 3 * 2, num_x // 3: num_x // 3 * 2]
        else:
            array_2d_filled = grid_z
    except ValueError:
        # The Richie Glitchie!
        msg = "  2D interpolation failed, most likely no good data"
        logger.warning(msg)
        array_2d_filled = array_2d.copy()
    # Check something...
    if WasMA:
        array_2d_filled = numpy.ma.masked_values(array_2d_filled, c.missing_value)
        array_2d_filled = numpy.ma.masked_invalid(array_2d_filled)
    # Return the filled array
    return array_2d_filled

def write_data_1columnpermonth(xlSheet, data, ts, format_string=''):
    xlCol = 0
    # write the data to the xl file
    nrows = numpy.shape(data)[0]
    ncols = numpy.shape(data)[1]
    xlSheet.write(1,xlCol,'Hour')
    for j in range(nrows+1):
        xlSheet.write(j+2,xlCol,float(j)*ts/60)
    xlCol = xlCol + 1
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,ncols+1):
        xlSheet.write(0,xlCol,calendar.month_abbr[m])
        xlSheet.write(1,xlCol,'Av')
        for j in range(nrows):
            xlSheet.write(j+2,xlCol,data[j,m-1],d_xf)
        xlCol = xlCol + 1

def write_data_1columnpertimestep(xlSheet, data, ts, startdate, format_string):
    tmp = data.copy()
    if numpy.ma.isMA(tmp): tmp = numpy.ma.filled(tmp,float(c.missing_value))
    xlCol = 0
    # write the data to the xl file
    xlSheet.write(1,xlCol,'Day')
    nrows = numpy.shape(tmp)[0]
    ncols = numpy.shape(tmp)[1]
    if startdate is None:
        for j in range(nrows+1):
            xlSheet.write(j+2,xlCol,j)
    else:
        d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy')
        for j in range(nrows):
            d = startdate + datetime.timedelta(days=j)
            xlSheet.write(j+2,xlCol,d,d_xf)
    xlCol = xlCol + 1
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,ncols+1):
        xlSheet.write(1,xlCol,float(m)*ts/60)
        for j in range(nrows):
            xlSheet.write(j+2,xlCol,float(tmp[j,m-1]),d_xf)
        xlCol = xlCol + 1

def do_diurnalstats(Month, Hdh, data, xlSheet, format_string='',ts=30):
    xlCol = 0
    nInts = 24*int((60/ts)+0.5)
    Av_all = numpy.ma.zeros([nInts,12]) + float(c.missing_value)
    if len(format_string)!=0:
        d_xf = xlwt.easyxf(num_format_str=format_string)
    else:
        d_xf = xlwt.easyxf()
    for m in range(1,13):
        mi = numpy.where(Month==m)[0]
        Num,Hr,Av,Sd,Mx,Mn = get_diurnalstats(Hdh[mi],data[mi],ts)
        Av_all[:,m-1] = Av[:]
        Num = numpy.ma.filled(Num,float(c.missing_value))
        Hr = numpy.ma.filled(Hr,float(c.missing_value))
        Av = numpy.ma.filled(Av,float(c.missing_value))
        Sd = numpy.ma.filled(Sd,float(c.missing_value))
        Mx = numpy.ma.filled(Mx,float(c.missing_value))
        Mn = numpy.ma.filled(Mn,float(c.missing_value))
        if m==1:
            xlSheet.write(1,xlCol,'Hour')
            for j in range(len(Hr)):
                xlSheet.write(j+2,xlCol,Hr[j])
            xlCol = xlCol + 1
        xlSheet.write(0,xlCol,calendar.month_abbr[m])
        xlSheet.write(1,xlCol,'Num')
        xlSheet.write(1,xlCol+1,'Av')
        xlSheet.write(1,xlCol+2,'Sd')
        xlSheet.write(1,xlCol+3,'Mx')
        xlSheet.write(1,xlCol+4,'Mn')
        for j in range(len(Hr)):
            xlSheet.write(j+2,xlCol,int(Num[j]))
            xlSheet.write(j+2,xlCol+1,Av[j],d_xf)
            xlSheet.write(j+2,xlCol+2,Sd[j],d_xf)
            xlSheet.write(j+2,xlCol+3,Mx[j],d_xf)
            xlSheet.write(j+2,xlCol+4,Mn[j],d_xf)
        xlCol = xlCol + 5
    return Av_all

def get_diurnalstats(DecHour,Data,ts):
    nInts = 24*int((60/ts)+0.5)
    Num = numpy.ma.zeros(nInts,dtype=int)
    Hr = numpy.ma.zeros(nInts,dtype=float)
    for i in range(nInts):
        Hr[i] = float(i)*ts/60.
    Av = numpy.ma.masked_all(nInts)
    Sd = numpy.ma.masked_all(nInts)
    Mx = numpy.ma.masked_all(nInts)
    Mn = numpy.ma.masked_all(nInts)
    if numpy.size(Data)!=0:
        for i in range(nInts):
            li = numpy.ma.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-float(c.missing_value))>c.eps))
            Num[i] = numpy.size(li)
            if Num[i]!=0:
                Av[i] = numpy.ma.mean(Data[li])
                Sd[i] = numpy.ma.std(Data[li])
                Mx[i] = numpy.ma.max(Data[li])
                Mn[i] = numpy.ma.min(Data[li])
    return Num, Hr, Av, Sd, Mx, Mn

def get_rangecheck_limit(cf,label,upr_def=1E10,lwr_def=-1E10):
    upper = float(upr_def)
    lower = float(lwr_def)
    for section in ['Variables']:
        if label in list(cf[section].keys()):
            if 'RangeCheck' in list(cf[section][label].keys()):
                upper = float(cf[section][label]['RangeCheck']['upper'])
                lower = float(cf[section][label]['RangeCheck']['lower'])
    return upper,lower

def get_formatstring(cf,label,fmt_def=''):
    fmt_str = fmt_def
    for section in ['Variables']:
        if label in list(cf[section].keys()):
            if 'format' in list(cf[section][label].keys()):
                fmt_str = str(cf[section][label]['format'])
    return fmt_str

def climatology(cf):
    nc_filename = pfp_io.get_infilenamefromcf(cf)
    if not pfp_utils.file_exists(nc_filename):
        msg = " Unable to find netCDF file " + nc_filename
        logger.error(msg)
        return
    xl_filename = nc_filename.replace(".nc", "_Climatology.xls")
    xlFile = xlwt.Workbook()
    ds = pfp_io.NetCDFRead(nc_filename)
    if ds.info["returncodes"]["value"] != 0:
        return
    # get the time step
    ts = int(ds.root["Attributes"]['time_step'])
    # get the datetime series
    dt = ds.root["Variables"]['DateTime']['Data']
    ldt = dt - datetime.timedelta(minutes=ts)
    start = datetime.datetime(ldt[0].year, ldt[0].month, ldt[0].day, 0, 0, 0)
    start += datetime.timedelta(minutes=ts)
    end = datetime.datetime(ldt[-1].year, ldt[-1].month, ldt[-1].day, 0, 0, 0)
    end += datetime.timedelta(minutes=1440)

    Hdh = numpy.array([(d.hour + d.minute/float(60)) for d in dt])
    Month = numpy.array([d.month for d in dt])
    # get the initial start and end dates
    StartDate = str(dt[0])
    EndDate = str(dt[-1])
    # find the start index of the first whole day (time=00:30)
    si = pfp_utils.GetDateIndex(dt,StartDate,ts=ts,default=0,match='startnextday')
    # find the end index of the last whole day (time=00:00)
    ei = pfp_utils.GetDateIndex(dt,EndDate,ts=ts,default=-1,match='endpreviousday')
    # get local views of the datetime series
    Hdh = Hdh[si:ei+1]
    Month = Month[si:ei+1]
    # get the number of time steps in a day and the number of days in the data
    ntsInDay = int(24.0*60.0/float(ts))
    # loop over the variables listed in the control file
    cf_labels = sorted(list(cf['Variables'].keys()))
    ds_labels = sorted(list(ds.root["Variables"].keys()))
    for label in cf_labels:
        # check to see if an alternative variable name is given
        label = pfp_utils.get_keyvaluefromcf(cf, ["Variables"], "name", default=label)
        if label in ds_labels:
            logger.info(" Doing climatology for " + label)
            var = pfp_utils.GetVariable(ds, label, start=si, end=ei)
            # do the diurnal by month statistics
            fmt_str = get_formatstring(cf, label, fmt_def='')
            xlSheet = xlFile.add_sheet(label)
            Av_all = do_diurnalstats(Month, Hdh, var["Data"], xlSheet, format_string=fmt_str, ts=ts)
            # do the daily statistics
            var = pfp_utils.GetVariable(ds, label)
            var = pfp_utils.PadVariable(var, start, end)
            nDays = int(len(var["Data"]))//ntsInDay
            data_daily = var["Data"].reshape(nDays, ntsInDay)
            xlSheet = xlFile.add_sheet(label + '(day)')
            write_data_1columnpertimestep(xlSheet, data_daily, ts, start, fmt_str)
            data_daily_i = do_2dinterpolation(data_daily)
            # check to see if the interpolation has left some data unfilled
            # this can happen if there is missing data at the boundaries of
            # the date x hour 2D array
            idx = numpy.where(numpy.ma.getmaskarray(data_daily_i) == True)
            # fill any missing data in the interpolated array with the monthly climatology
            month = numpy.array([dt.month-1 for dt in var["DateTime"]])
            month_daily = month.reshape(nDays, ntsInDay)
            hour = numpy.array([int(60/ts*(dt.hour+float(dt.minute)/float(60)))
                                for dt in var["DateTime"]])
            hour_daily = hour.reshape(nDays, ntsInDay)
            data_daily_i[idx] = Av_all[hour_daily[idx], month_daily[idx]]
            # write the interpolated data to the Excel workbook
            xlSheet = xlFile.add_sheet(label + 'i(day)')
            write_data_1columnpertimestep(xlSheet, data_daily_i, ts, start, fmt_str)
        else:
            msg = " Requested variable " + label + " not in data structure"
            logger.warning(msg)
            continue
    msg = " Saving Excel file " + os.path.split(xl_filename)[1]
    logger.info(msg)
    xlFile.save(xl_filename)
    return
