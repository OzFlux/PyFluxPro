# standard modules
import copy
import datetime
import faulthandler
import logging
import os
import shutil
import sys
import _pickle as cPickle
# 3rd party modules
from configobj import ConfigObj
import cftime
import dateutil
import netCDF4
import numpy
import xlrd
# PFP modules
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("external")-1])
import scripts.constants as c
import scripts.meteorologicalfunctions as pfp_mf
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log
import scripts.pfp_ts as pfp_ts
import scripts.pfp_utils as pfp_utils

faulthandler.enable()

now = datetime.datetime.now()
log_file_name = "access2nc_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_path = "logfiles"
if not os.path.isdir(log_file_path):
    os.mkdir(log_file_path)
log_file_uri = os.path.join(log_file_path, log_file_name)
logger = pfp_log.CreateLogger("pfp_log", log_file_name=log_file_uri, to_screen=True)

def read_site_master(site_master_file_path, sheet_name):
    """
    """
    xl_book = xlrd.open_workbook(site_master_file_path)
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
    site_info = {}
    for n in range(first_data_row,last_row):
        site_name = xl_sheet.cell(n,0).value
        site_name = site_name.replace(" ","")
        site_info[site_name] = {}
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value
    return site_info
def build_concatenation_dictionary(site_info, cfg):
    existing_access_base_path = cfg["Files"]["existing_access_base_path"]
    concatenation_info = {}
    sites = list(site_info.keys())
    for site in sites:
        access_file_path = os.path.join(existing_access_base_path, site, "Data", "ACCESS")
        access_file_name = site + "_ACCESS" + ".nc"
        out_file_name = os.path.join(access_file_path, access_file_name)
        if cfg["Files"].as_bool("from_cleaned"):
            access_file_path = os.path.join(access_file_path, "cleaned")
        in_file_name = os.path.join(access_file_path, access_file_name)
        in_file_names = [in_file_name]
        concatenation_info[site] = {"NetCDFConcatenate": {}, "RemoveIntermediateSeries": {}}
        cisn = concatenation_info[site]["NetCDFConcatenate"]
        cisn["OK"] = True
        cisn["in_file_names"] = in_file_names
        cisn["out_file_name"] = out_file_name
        cisn["NumberOfDimensions"] = 3
        cisn["MaxGapInterpolate"] = 0
        cisn["FixTimeStepMethod"] = "round"
        cisn["Truncate"] = "No"
        cisn["TruncateThreshold"] = 50
        cisn["SeriesToCheck"] = ["all"]
        cisn["time_coverage_start"] = []
        cisn["time_coverage_end"] = []
        cisn["chrono_files"] = []
        cisn["labels"] = []
        cisn["attributes"] = ['long_name', 'standard_name', 'statistic_type', 'units']
        cisr = concatenation_info[site]["RemoveIntermediateSeries"]
        cisr["KeepIntermediateSeries"] = "No"
        cisr["not_output"] = []
    return concatenation_info
def convert_utc_to_local_standard(utc_notz, site_timezone):
    # get the UTC timezone object
    from_zone = dateutil.tz.gettz('UTC')
    # get the site timezone object
    to_zone = dateutil.tz.gettz(site_timezone)
    # add the UTC timezone information to the ACCESS datetime
    utc_tz = [dt.replace(tzinfo=from_zone) for dt in utc_notz]
    # get the local daylight savings time at the site from UTC
    loc_dst = [dt.astimezone(to_zone) for dt in utc_tz]
    # convert from daylight savings time to standard time at the site
    loc_ast = [dt-dt.dst() for dt in loc_dst]
    # remove the timezone information from the local standard time at the site
    loc_ast_notz = [dt.replace(tzinfo=None) for dt in loc_ast]
    return loc_ast_notz
def delete_intermediate_variables(ds, intermediate_variables):
    labels = list(ds.series.keys())
    for iv_label in intermediate_variables:
        for i in range(3):
            for j in range(3):
                label = iv_label+"_"+str(i)+str(j)
                if label in labels:
                    ds.series.pop(label)
    return
def init_data(cfg, site_info):
    labels = list(cfg["Variables"].keys())
    sites = list(site_info.keys())
    # initialise the data dictionary
    data = {}
    for site in sites:
        data[site] = {}
        data[site]["file_id"] = []
        data[site]["dt_utc"] = []
        data[site]["DateTime"] = {"Data": [],
                                  "Flag": [],
                                  "Attr": {"long_name": "Datetime in local timezone",
                                           "units": ""}}
        for label in labels:
            attr = copy.deepcopy(cfg["Variables"][label]["Attr"])
            for i in range(3):
                for j in range(3):
                    # we use lists because they are fast
                    data[site][label+"_"+str(i)+str(j)] = {"Data": [],
                                                           "Flag": [],
                                                           "Attr": attr}
    return data
def calculate_upwelling_shortwave_radiation(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            Fsd = pfp_utils.GetVariable(ds, "Fsd"+"_"+str(i)+str(j))
            Fnsw = pfp_utils.GetVariable(ds, "Fnsw"+"_"+str(i)+str(j))
            Fsu = pfp_utils.CreateEmptyVariable("Fsu"+"_"+str(i)+str(j), nrecs,
                                                datetime=Fsd["DateTime"])
            Fsu["Data"] = Fsd["Data"] - Fnsw["Data"]
            Fsu["Flag"] = numpy.where(numpy.ma.getmaskarray(Fsd["Data"]) == True, ones, zeros)
            Fsu["Attr"] = {"standard_name": "surface_upwelling_shortwave_flux_in_air",
                           "long_name": "Up-welling shortwave radiation",
                           "units": "W/m^2", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Fsu)
    return
def calculate_upwelling_longwave_radiation(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            Fld = pfp_utils.GetVariable(ds, "Fld"+"_"+str(i)+str(j))
            Fnlw = pfp_utils.GetVariable(ds, "Fnlw"+"_"+str(i)+str(j))
            Flu = pfp_utils.CreateEmptyVariable("Flu"+"_"+str(i)+str(j), nrecs,
                                                datetime=Fld["DateTime"])
            Flu["Data"] = Fld["Data"] - Fnlw["Data"]
            Flu["Flag"] = numpy.where(numpy.ma.getmaskarray(Fld["Data"]) == True, ones, zeros)
            Flu["Attr"] = {"standard_name": "surface_upwelling_longwave_flux_in_air",
                           "long_name": "Up-welling shortwave radiation",
                           "units": "W/m^2", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Flu)
    return
def calculate_net_radiation(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            Fnsw = pfp_utils.GetVariable(ds, "Fnsw"+"_"+str(i)+str(j))
            Fnlw = pfp_utils.GetVariable(ds, "Fnlw"+"_"+str(i)+str(j))
            Fn = pfp_utils.CreateEmptyVariable("Fn"+"_"+str(i)+str(j), nrecs,
                                                datetime=Fnsw["DateTime"])
            Fn["Data"] = Fnsw["Data"] + Fnlw["Data"]
            Fn["Flag"] = numpy.where(numpy.ma.getmaskarray(Fn["Data"]) == True, ones, zeros)
            Fn["Attr"] = {"standard_name": "surface_net_downwawrd_radiative_flux",
                           "long_name": "Net radiation",
                           "units": "W/m^2", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Fn)
    return
def calculate_available_energy(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            Fh = pfp_utils.GetVariable(ds, "Fh"+"_"+str(i)+str(j))
            Fe = pfp_utils.GetVariable(ds, "Fe"+"_"+str(i)+str(j))
            Fa = pfp_utils.CreateEmptyVariable("Fa"+"_"+str(i)+str(j), nrecs,
                                                datetime=Fh["DateTime"])
            Fa["Data"] = Fh["Data"] + Fe["Data"]
            Fa["Flag"] = numpy.where(numpy.ma.getmaskarray(Fa["Data"]) == True, ones, zeros)
            Fa["Attr"] = {"long_name": "Available energy",
                          "units": "W/m^2", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Fa)
    return
def calculate_ground_heat_flux(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            Fn = pfp_utils.GetVariable(ds, "Fn"+"_"+str(i)+str(j))
            Fa = pfp_utils.GetVariable(ds, "Fa"+"_"+str(i)+str(j))
            Fg = pfp_utils.CreateEmptyVariable("Fg"+"_"+str(i)+str(j), nrecs,
                                                datetime=Fn["DateTime"])
            Fg["Data"] = Fn["Data"] - Fa["Data"]
            Fg["Flag"] = numpy.where(numpy.ma.getmaskarray(Fg["Data"]) == True, ones, zeros)
            Fg["Attr"] = {"standard_name": "downward_heat_flux_in_soil",
                           "long_name": "Ground heat flux",
                           "units": "W/m^2", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Fg)
    return
def calculate_absolute_humidity(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            RH = pfp_utils.GetVariable(ds, "RH"+"_"+str(i)+str(j))
            Ta = pfp_utils.GetVariable(ds, "Ta"+"_"+str(i)+str(j))
            AH = pfp_utils.CreateEmptyVariable("AH"+"_"+str(i)+str(j), nrecs,
                                                datetime=RH["DateTime"])
            AH["Data"] = pfp_mf.absolutehumidityfromrelativehumidity(Ta["Data"], RH["Data"])
            AH["Flag"] = numpy.where(numpy.ma.getmaskarray(AH["Data"]) == True, ones, zeros)
            AH["Attr"] = {"standard_name": "mass_concentration_of_water_vapor_in_air",
                           "long_name": "Absolute humidity",
                           "units": "g/m^3", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, AH)
    return
def calculate_ws_and_wd_from_u_and_v(ds):
    nrecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    for i in range(3):
        for j in range(3):
            U = pfp_utils.GetVariable(ds, "U"+"_"+str(i)+str(j))
            V = pfp_utils.GetVariable(ds, "V"+"_"+str(i)+str(j))
            Ws = pfp_utils.CreateEmptyVariable("Ws"+"_"+str(i)+str(j), nrecs,
                                                datetime=U["DateTime"])
            Wd = pfp_utils.CreateEmptyVariable("Wd"+"_"+str(i)+str(j), nrecs,
                                                datetime=U["DateTime"])
            # get the wind speed and direction from the components
            Wd["Data"] = float(270) - (numpy.degrees(numpy.ma.arctan2(V["Data"], U["Data"])))
            Wd["Data"] = numpy.ma.mod(Wd["Data"], 360)
            Ws["Data"] = numpy.ma.sqrt(U["Data"]*U["Data"] + V["Data"]*V["Data"])
            # mask wind direction when the wind speed is less than 0.01
            Wd["Data"] = numpy.ma.masked_where(Ws["Data"] < 0.01, Wd["Data"])
            # now set the QC flag
            Ws["Flag"] = numpy.where(numpy.ma.getmaskarray(Ws["Data"]) == True, ones, zeros)
            Wd["Flag"] = numpy.where(numpy.ma.getmaskarray(Wd["Data"]) == True, ones, zeros)
            # update the variable attributes
            Ws["Attr"] = {"standard_name": "wind_speed", "long_name": "Wind speed",
                          "units": "m/s", "statistic_type": "average"}
            Wd["Attr"] = {"standard_name": "wind_from_direction", "long_name": "Wind direction",
                          "units": "degrees", "statistic_type": "average"}
            pfp_utils.CreateVariable(ds, Ws)
            pfp_utils.CreateVariable(ds, Wd)
    return

logger = logging.getLogger("pfp_log")
# read the control file
cfg_file_path = "process_access2nc.txt"
msg = " Loading the control file"
logger.info(msg)
cfg = pfp_io.get_controlfilecontents(cfg_file_path)
cfg_labels = [l for l in list(cfg["Variables"].keys()) if "nc" in list(cfg["Variables"][l].keys())]
# read the site master workbook
site_info = read_site_master(cfg["Files"]["site_master_file_path"], "ACCESS")
sites = list(site_info.keys())
#sites = ["Calperum"]
# build the information dictionary for concatenation
concatenation_info = build_concatenation_dictionary(site_info, cfg)

dt_utc = []
data = init_data(cfg, site_info)
new_access_base_path = cfg["Files"]["new_access_base_path"]
new_access_directories = sorted([f.path for f in os.scandir(new_access_base_path) if f.is_dir()])
n_new_access_directories = len(new_access_directories)
if n_new_access_directories == 0:
    msg = " No files found in " + new_access_base_path
    logger.info(msg)
    sys.exit()

if cfg["Files"].as_bool("from_file"):
    # loop over the ACCESS directories, read the ACCESS files, extract a 3 x 3 cut out around the site
    # and put the data into a dictionary ordered by site
    msg = " Processing ACCESS files from " + new_access_directories[0].split(os.sep)[-1]
    msg += " to " + new_access_directories[-1].split(os.sep)[-1]
    logger.info(msg)
    warnings = {"all": [], "cant_read": [], "netcdf_error": []}
    for new_access_directory in new_access_directories:
        fs = os.scandir(new_access_directory)
        new_access_files = sorted([f.name for f in fs if "IDY25026" in f.name])
        n_new_access_files = len(new_access_files)
        for new_access_file in new_access_files:
            # get the file ID, this is YYYYMMDDHHnnn where
            #  HH is the base time e.g. 00, 06, 12, 18
            #  nnn is the forecast time e.g. 000, 001, 002, 003, 004, 005, 006
            parts = new_access_file.split(".")
            file_id = parts[4] + "_" + parts[5]
            new_access_file_path = os.path.join(new_access_directory, new_access_file)
            try:
                access_file = netCDF4.Dataset(new_access_file_path)
            except:
                base_name = os.path.basename(new_access_file_path)
                msg = " Can't read " + base_name + ", skipping file ..."
                #logger.warning(msg)
                warnings["cant_read"].append(msg)
                #access_file.close()
                continue
            # get the Python datetime (UTC)
            nc_time_data = numpy.ma.filled(access_file.variables["time"])[0]
            nc_time_units = getattr(access_file.variables["time"], "units")
            dt_utc = cftime.num2pydate(nc_time_data, nc_time_units)
            # get the latitude and longitude
            latitude = access_file.variables["lat"]
            longitude = access_file.variables["lon"]
            # get the latitude and longitude resolution
            delta_latitude = (max(latitude) - min(latitude))/(len(latitude)-1)
            delta_longitude = (max(longitude) - min(longitude))/(len(longitude)-1)
            # loop over sites
            for site in sites:
                # get the site latitude and longitude
                site_latitude = site_info[site]["Latitude"]
                site_longitude = site_info[site]["Longitude"]
                #
                # we need a check here to make sure the selected grid square is land, not water
                #
                y = int(((latitude[0]-site_latitude)/delta_latitude)+0.5)
                x = int(((site_longitude-longitude[0])/delta_longitude)+0.5)
                # loop over variables to read from the ACCESS file
                data[site]["file_id"].append(file_id)
                data[site]["dt_utc"].append(dt_utc)
                n_cfg_labels = len(cfg_labels)
                for label in cfg_labels:
                    d = numpy.full([3, 3], -9999, dtype=numpy.float32)
                    access_label = cfg["Variables"][label]["nc"]["name"]
                    access_variable = access_file.variables[access_label]
                    # check the variable dimensions
                    if len(access_variable.shape) == 3:
                        # most ACCESS variables have 3 dimensions
                        try:
                            d = access_variable[0, y-1:y+2, x-1:x+2].data
                        except:
                            base_name = os.path.basename(new_access_file_path)
                            msg = " netCDF error: " + base_name + " (" + site + ", "
                            msg += label + ")"
                            #logger.warning(msg)
                            warnings["netcdf_error"].append(msg)
                            # create dummy variable values
                            d = numpy.full((3, 3), float(c.missing_value), dtype=numpy.float32)
                    elif len(access_variable.shape) == 4:
                        # but soil moisture and temperature have 4 ...
                        try:
                            d = access_variable[0, 0, y-1:y+2, x-1:x+2].data
                        except:
                            base_name = os.path.basename(new_access_file_path)
                            msg = " netCDF error: " + base_name + " (" + site + ", "
                            msg += label + ")"
                            #logger.warning(msg)
                            warnings["netcdf_error"].append(msg)
                            # create dummy variable values
                            d = numpy.full((3, 3), float(c.missing_value), dtype=numpy.float32)
                    else:
                        msg = " Unexpected dimensions in netCDF file"
                        logger.error(msg)
                        raise Exception
                    # loop over the 3 x 3 cut out around the site
                    for i in range(3):
                        for j in range(3):
                            data[site][label+"_"+str(i)+str(j)]["Data"].append(d[i, j])
            access_file.close()
        # move the processed files from "download" to "processed"
        msg = " Moving " + new_access_directory.split(os.sep)[-1] + " to processed"
        logger.info(msg)
        #processed_access_directory = new_access_directory.replace("downloaded", "processed")
        #shutil.move(new_access_directory, processed_access_directory)
    wlist = [warnings[k] for k in warnings.keys()]
    warnings["all"] = [item for sublist in wlist for item in sublist]
    if len(warnings["all"]) != 0:
        n_warnings = len(list(set(warnings["all"])))
        msg = str(n_warnings) + " warnings occurred"
        logger.warning("")
        logger.warning(msg)
        logger.warning("")
    # pickle the data in case we want to use it again
    with open(r"data.pickle", "wb") as output_file:
        cPickle.dump(data, output_file)
else:
    # read the pickled data
    msg = " Reading from pickle file"
    logger.info(msg)
    import _pickle as cPickle
    with open(r"data.pickle", "rb") as input_file:
        data = cPickle.load(input_file)

# create a dictionary of data sets, one for each site, with no gaps and copy the data into it
# first we create a list of all posible file IDs between the first
# one read and the last one read
date_string = os.path.normpath(new_access_directories[0]).split(os.path.sep)[-1]
start = datetime.datetime.strptime(date_string, "%Y%m%d%H")
date_string = os.path.normpath(new_access_directories[-1]).split(os.path.sep)[-1]
end = datetime.datetime.strptime(date_string, "%Y%m%d%H")
delta_dt = datetime.timedelta(hours=6)
dt_6hourly = [dt for dt in pfp_utils.perdelta(start, end, delta_dt)]
file_id_nogaps = []
dt_utc_nogaps = []
for dt in dt_6hourly:
    for n in ["_000", "_001", "_002", "_003", "_004", "_005", "_006"]:
        dts = dt.strftime("%Y%m%d%H") + n
        file_id_nogaps.append(dts)
        parts = dts.split("_")
        dtd = datetime.datetime.strptime(parts[0], "%Y%m%d%H") + datetime.timedelta(hours=int(parts[1]))
        dt_utc_nogaps.append(dtd)
dt_utc_nogaps = numpy.array(dt_utc_nogaps)
# number of records in the no-gap data set
nrecs_nogaps = len(file_id_nogaps)
# check to see if both lists contain the same elements
# create a data set for the no gap data
data_nogaps = init_data(cfg, site_info)
msg = " Creating data sets without gaps"
logger.info(msg)
for site in sites:
    # put the file_id and the UTC datetime into the no-gaps data structure
    data_nogaps[site]["file_id"] = file_id_nogaps
    data_nogaps[site]["dt_utc"] = dt_utc_nogaps

    # get the set of elements common to both file_id_nogaps and the file_id for this site
    both = set(file_id_nogaps).intersection(data[site]["file_id"])
    # get the indices of common elements in file_id_nogaps
    iA = [file_id_nogaps.index(x) for x in both]
    # get the indices of common elements in the site file_id
    iB = [data[site]["file_id"].index(x) for x in both]

    # do the datetime
    data_nogaps[site]["DateTime"]["Data"] = dt_utc_nogaps
    data_nogaps[site]["DateTime"]["Flag"] = numpy.zeros(nrecs_nogaps)
    data_nogaps[site]["DateTime"]["Attr"] = {"long_name": "Datetime in UTC timezone", "units": ""}

    # loop over variables to read from the ACCESS file
    for label in cfg_labels:
        for i in range(3):
            for j in range(3):
                l = label+"_"+str(i)+str(j)
                data_nogaps[site][l]["Data"] = numpy.ma.masked_all(nrecs_nogaps)
                data_nogaps[site][l]["Flag"] = numpy.ones(nrecs_nogaps)
                data_nogaps[site][l]["Data"][iA] = numpy.ma.array(data[site][l]["Data"])[iB]
                data_nogaps[site][l]["Flag"][iA] = int(0)

# now we copy the data from the no gaps data sets to PFP data structures
# dictionary to hold data structures for each site
dss_ats = {}
msg = " Creating data structure for each site"
logger.info(msg)
for site in sites:
    # create a data structure for this site
    dss_ats[site] = pfp_io.DataStructure()
    # convert UTC datetime to local standard time
    dt_loc_nogaps = numpy.array(convert_utc_to_local_standard(dt_utc_nogaps, site_info[site]["Time zone"]))
    # add the global attributes
    dss_ats[site].globalattributes["site_name"] = site.replace(" ", "")
    dss_ats[site].globalattributes["time_zone"] = site_info[site]["Time zone"]
    dss_ats[site].globalattributes["latitude"] = site_info[site]["Latitude"]
    dss_ats[site].globalattributes["longitude"] = site_info[site]["Longitude"]
    dss_ats[site].globalattributes["time_coverage_start"] = str(dt_loc_nogaps[0])
    dss_ats[site].globalattributes["time_coverage_end"] = str(dt_loc_nogaps[-1])
    dss_ats[site].globalattributes["processing_level"] = "L1"
    # we discard the base time plus 6 hour data for non-precipitation variables
    file_id = data_nogaps[site]["file_id"]
    idx = numpy.array([i for i, l in enumerate(file_id) if l[-3:] != "006"])
    nrecs = len(idx)
    dss_ats[site].globalattributes["nc_nrecs"] = nrecs
    zeros = numpy.zeros(nrecs)
    ones = numpy.ones(nrecs)
    # put the Python datetime in the data structure
    attr = {"long_name": "Datetime in local timezone", "units": ""}
    variable = {"Label": "DateTime", "Data": dt_loc_nogaps[idx], "Flag": zeros, "Attr": attr}
    pfp_utils.CreateVariable(dss_ats[site], variable)
    # get the time step
    ts = numpy.mean(pfp_utils.get_timestep(dss_ats[site])//60)
    ts = pfp_utils.roundtobase(ts, base=30)
    dss_ats[site].globalattributes["time_step"] = str(int(ts))
    # put the UTC datetime in the data structure
    attr = {"long_name": "Datetime in UTC", "units": ""}
    variable = {"Label": "DateTime_UTC", "Data": dt_utc_nogaps[idx], "Flag": zeros, "Attr": attr}
    pfp_utils.CreateVariable(dss_ats[site], variable)
    # loop over the non-precipitation variables
    these_labels = [l for l in cfg_labels if "Precip" not in l]
    for label in these_labels:
        for i in range(3):
            for j in range(3):
                d = data_nogaps[site][label+"_"+str(i)+str(j)]["Data"][idx]
                # create the variable
                variable = {"Label": label+"_"+str(i)+str(j),
                            "Data": d,
                            "Flag": numpy.where(numpy.ma.getmaskarray(d) == True, ones, zeros),
                            "Attr": data[site][label+"_"+str(i)+str(j)]["Attr"]}
                # put it into the data structure
                pfp_utils.CreateVariable(dss_ats[site], variable)

# do the units conversion
msg = " Converting units for ps, Ta, Ts and Sws"
logger.info(msg)
for site in sites:
    for i in range(3):
        for j in range(3):
            # surface pressure units Pa to kPa
            ps = pfp_utils.GetVariable(dss_ats[site], "ps"+"_"+str(i)+str(j))
            ps = pfp_utils.convert_units_func(dss_ats[site], ps, "kPa")
            pfp_utils.CreateVariable(dss_ats[site], ps)
            # air temperature
            Ta = pfp_utils.GetVariable(dss_ats[site], "Ta"+"_"+str(i)+str(j))
            Ta = pfp_utils.convert_units_func(dss_ats[site], Ta, "degC")
            pfp_utils.CreateVariable(dss_ats[site], Ta)
            # soil temperature
            Ts = pfp_utils.GetVariable(dss_ats[site], "Ts"+"_"+str(i)+str(j))
            Ts = pfp_utils.convert_units_func(dss_ats[site], Ts, "degC")
            pfp_utils.CreateVariable(dss_ats[site], Ts)
            # soil moisture
            Sws = pfp_utils.GetVariable(dss_ats[site], "Sws"+"_"+str(i)+str(j))
            Sws = pfp_utils.convert_units_func(dss_ats[site], Sws, "m^3/m^3")
            pfp_utils.CreateVariable(dss_ats[site], Sws)

# get the derived quantities
msg = " Deriving radiation, ground heat flux, absolute humidity, wind speed and direction"
logger.info(msg)
for site in sites:
    calculate_upwelling_shortwave_radiation(dss_ats[site])
    calculate_upwelling_longwave_radiation(dss_ats[site])
    calculate_net_radiation(dss_ats[site])
    calculate_available_energy(dss_ats[site])
    calculate_ground_heat_flux(dss_ats[site])
    calculate_absolute_humidity(dss_ats[site])
    calculate_ws_and_wd_from_u_and_v(dss_ats[site])
    delete_intermediate_variables(dss_ats[site], ["Fnlw", "Fnsw"])

# now do the precipitation
msg = " Doing precipitation each site "
logger.info(msg)
for site in sites:
    # precipitation from ACCESS is an accumulated quantity which means that the instantaneous
    # precip at the model base times (00, 06, 12, 18) can't be calculated (we don't have data
    # for base time minus 1 hour).  So, we use the difference in accumulated precip between
    # the 5th and 6th forecast hour from the previous base time as the instantaneous precip
    # for the following base time.
    those_labels = [l for l in cfg_labels if "Precip" in l]
    file_id = data_nogaps[site]["file_id"]
    idx = numpy.array([i for i, l in enumerate(file_id) if l[-3:] != "006"])
    for label in those_labels:
        for i in range(3):
            for j in range(3):
                # accummulated precipitation
                ap = data_nogaps[site][label+"_"+str(i)+str(j)]["Data"]
                # instantaneous precipitation
                ip = numpy.ma.ediff1d(ap, to_begin=ap[0])
                # shift instantaneous precipitation 1 place to the right
                sip = ip.copy()
                sip[0] = ip[-1]
                sip[1:] = ip[:-1]
                # get the indices of files ending in "000"
                idx_000 = numpy.array([i for i, l in enumerate(file_id) if l[-3:] == "000"])
                # put the 006-005 precipitation into 000
                ip[idx_000] = sip[idx_000]
                d = ip[idx]
                # create the variable
                attr = data[site][label+"_"+str(i)+str(j)]["Attr"]
                attr["statistic_type"] = "sum"
                variable = {"Label": label+"_"+str(i)+str(j),
                            "Data": d,
                            "Flag": numpy.where(numpy.ma.getmaskarray(d) == True, ones, zeros),
                            "Attr": attr}
                # put it into the data structure
                pfp_utils.CreateVariable(dss_ats[site], variable)

# interpolate from the ACCESS time step (60 minutes) to the tower time step as required
dss_tts = {}
msg = " Interpolating data for each site "
logger.info(msg)
for site in sites:
    access_ts = int(float(dss_ats[site].globalattributes["time_step"]))
    tower_ts = int(float(site_info[site]["Time step"]))
    labels = [l for l in list(dss_ats[site].series.keys()) if "DateTime" not in l]
    # interpolate from the ACCESS time step (60 minutes) to the tower time step
    dss_tts[site] = pfp_ts.InterpolateDataStructure(dss_ats[site], labels=labels,
                                                    new_time_step=tower_ts,
                                                    sums="interpolate",
                                                    mode="quiet")

# write out the files of ACCESS data
existing_access_base_path = cfg["Files"]["existing_access_base_path"]
msg = " Writing ACCESS netCDF files"
logger.info(msg)
for site in sites:
    cis = concatenation_info[site]
    # get the ACCESS data path
    access_file_path = os.path.join(existing_access_base_path, site, "Data", "ACCESS", "processed")
    # check the directory exists and create if it doesn't
    if not os.path.isdir(access_file_path):
        os.makedirs(access_file_path)
    # build the ACCESS file name
    dt = pfp_utils.GetVariable(dss_tts[site], "DateTime")
    start_date = dt["Data"][0].strftime("%Y%m%d%H%M")
    end_date = dt["Data"][-1].strftime("%Y%m%d%H%M")
    access_file_name = site + "_ACCESS" + "_" + start_date + "_" + end_date + ".nc"
    # get the full path including the file name
    access_file_uri = os.path.join(access_file_path, access_file_name)
    # and write the ACCESS data to a netCDF file
    pfp_io.NetCDFWrite(access_file_uri, dss_tts[site])
    cis["NetCDFConcatenate"]["in_file_names"].append(access_file_uri)
    # concatenate with the existing ACCESS data
    pfp_io.NetCDFConcatenate(cis)
    # save the concatenation_info dictionary
    cis_file_name = "cis_" + site + ".txt"
    cis_file_path = os.path.join(existing_access_base_path, site, "Data", "ACCESS", "cis")
    if not os.path.isdir(cis_file_path):
        os.mkdir(cis_file_path)
    cfg_cis = ConfigObj(cis, indent_type="    ", list_values=False)
    cfg_cis.filename = os.path.join(cis_file_path, cis_file_name)
    cfg_cis.write()

logger.info("All finished")