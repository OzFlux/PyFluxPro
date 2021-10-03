#!/usr/bin/env python3
"""
 Purpose:
  Download requested ISD files from the NOAA FTP site.
 Usage:
  python3 get_isd_from_noaa.py <controlfile>
 Author: PRI
 Date: June 2017
 modified December 2020 CEwenz
"""
# standard
from collections import OrderedDict
import datetime
import ftplib
import logging
import os
import sys
import time
# 3rd party
from configobj import ConfigObj
import xlrd
# check the scripts directory is present
#if not os.path.exists("scripts"):
#    print("get_isd_from_noaa: the scripts directory is missing")
#    sys.exit()
# since the scripts directory is there, try importing the modules
#sys.path.append('scripts')
# PFP
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("external")-1])
import scripts.pfp_io as pfp_io
import scripts.pfp_log as pfp_log

dir_list = ["logfiles/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)

now = datetime.datetime.now()
log_file_name = 'get_isd_from_noaa_' + now.strftime("%Y%m%d%H%M") + ".log"
log_file_name = os.path.join("logfiles", log_file_name)
logger = pfp_log.init_logger("pfp_log", log_file_name, to_file=False, to_screen=True)

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
        for item in header_row_values[1:]:
            i = header_row_values.index(item)
            site_info[site_name][item] = xl_sheet.cell(n,i).value

    return site_info

# read the control file file
if (__name__ == '__main__'):
    # get the control file name
    if len(sys.argv) == 1:
        # not on the command line, so ask the user
        cfg_file_path = input("Enter the control file name: ")
        # exit if nothing selected
        if len(cfg_file_path) == 0:
            sys.exit()
    else:
        # control file name on the command line
        if not os.path.exists(sys.argv[1]):
            # control file doesn't exist
            logger.error("Control file %s does not exist", sys.argv[1])
            sys.exit()
        else:
            cfg_file_path = sys.argv[1]

cf = pfp_io.get_controlfilecontents(cfg_file_path, mode="verbose")
xl_file_path  = cf["Files"]["xl_file_path"]
xl_sheet_name = cf["Files"]["xl_sheet_name"]
isd_base_path = cf["Files"]["isd_base_path"]
# get the site information from the site master spreadsheet
logger.info("Processing site master spreadsheet")
site_info = read_site_master(xl_file_path, xl_sheet_name)
# get a list of sites
site_list = list(site_info.keys())
# parse the control file to get a dictionary containing lists of the ISD sites
# required for each year
ftp_site_list = {}
for site in site_list:
    # get the list of ISD stations to be used for this site
    isd_site_list = []
    for item in ["ISD_ID_1","ISD_ID_2","ISD_ID_3","ISD_ID_4"]:
        if len(site_info[site][item]) != 0:
            isd_site_list.append(site_info[site][item])
    start_year = int(site_info[site]["Start year"])
    try:
        end_year = int(site_info[site]["End year"])
    except:
        end_year = int(now.year)
    year_list = [str(yr) for yr in range(start_year,end_year+1)]
    for year in year_list:
        if year not in ftp_site_list:
            ftp_site_list[year] = []
        ftp_site_list[year] += isd_site_list
# remove duplicate ISD site entries
year_list = sorted(ftp_site_list.keys())
for year in year_list:
    ftp_site_list[year] = list(set(ftp_site_list[year]))
# connect to the NOAA ISD ftp site
logger.info("Opening NOAA FTP site")
ftp_base_path = os.path.join("pub","data","noaa")
ftp = ftplib.FTP("ftp.ncdc.noaa.gov","ftp","pisaac.ozflux@gmail.com")
# loop over the years and get the ISD files
for year in year_list:
    logger.info("Processing year "+year)
    # get the directory for this year on the FTP server
    ftp_year_dir = os.path.join(ftp_base_path,year)
    # get a list of site files in this directory
    logger.info("Getting the ISD file list for "+year)
    ftp_file_list = ftp.nlst(ftp_year_dir)
    # get the path where site files from this year will be written
    out_year_path = os.path.join(isd_base_path,year)
    # create the directory if it doesn't exist
    if not os.path.exists(out_year_path):
        os.makedirs(out_year_path)
    # now loop over the ISD sites required for this year
    for isd_site in ftp_site_list[year]:
        # construct the ISD file name
        ftp_file_name = isd_site+"-"+str(year)+".gz"
        # and get the full path
        ftp_file_path = os.path.join(ftp_year_dir,ftp_file_name)
        # check that the site file exists in the FTP directory
        if ftp_file_path in ftp_file_list:
            # if it does, download it
            logger.info("Getting data for "+isd_site+" for year "+year)
            out_file_path = os.path.join(out_year_path,ftp_file_name)
            logger.info("Writing to "+out_file_path)
            ftp.retrbinary("RETR " + ftp_file_path, open(out_file_path, 'wb').write)
        else:
            logger.info("ISD site "+isd_site+" not found for year "+year)
# say goodbye
ftp.quit()
