# standard modules
import datetime
import faulthandler
import logging
import os
import platform
import sys
import traceback
import warnings
# 3rd party modules
from configobj import ConfigObj
import matplotlib
import netCDF4
# check the OS and set the matplotlib backend accordingly
if platform.system() == "Darwin":
    # set backend to "macosx" on Macs
    matplotlib.use("macosx")
elif platform.system() == "Windows":
    # set backend to "QT5Agg" for Windows
    matplotlib.use("QT5Agg")
elif platform.system() == "Linux":
    # set backend to "QT5Agg" for Linux
    matplotlib.use("QT5Agg")
else:
    # use whatever ...
    pass
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
# PFP modules
sys.path.append('scripts')
from scripts import cfg
from scripts import pfp_compliance
from scripts import pfp_gui
from scripts import pfp_io
from scripts import pfp_log
from scripts import pfp_threading
from scripts import pfp_top_level

faulthandler.enable()
warnings.filterwarnings("ignore", category=Warning)

# now check the logfiles and plots directories are present
dir_list = ["logfiles/", "plots/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# now check the solo/inf, solo/input, solo/log and solo/output directories are present
dir_list = ["solo/inf/", "solo/input/", "solo/log/", "solo/output/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# next we make sure the MPT directories are present ...
dir_list = ["mpt/input/", "mpt/log/", "mpt/output/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)
# ... and make sure the MDS directories are present
dir_list = ["mds/input/", "mds/log/", "mds/output/"]
for item in dir_list:
    if not os.path.exists(item):
        os.makedirs(item)

now = datetime.datetime.now()
logger_name = "pfp_log"
log_file_name = "pfp_" + now.strftime("%Y%m%d%H%M") + ".log"
log_file_name = os.path.join("logfiles", log_file_name)
logger = pfp_log.CreateLogger(logger_name, log_file_name=log_file_name)

class pfp_main_ui(QWidget):
    def __init__(self, pfp_version, textBox):
        super(pfp_main_ui, self).__init__()
        # menu bar
        self.menubar = QMenuBar(self)
        # File menu
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setTitle("File")
        # File/Convert submenu
        self.menuFileConvert = QMenu(self.menuFile)
        self.menuFileConvert.setTitle("Convert")
        # Edit menu
        self.menuEdit = QMenu(self.menubar)
        self.menuEdit.setTitle("Edit")
        # Run menu
        self.menuRun = QMenu(self.menubar)
        self.menuRun.setTitle("Run")
        # Plot menu
        self.menuPlot = QMenu(self.menubar)
        self.menuPlot.setTitle("Plot")
        # Utilities menu
        self.menuUtilities = QMenu(self.menubar)
        self.menuUtilities.setTitle("Utilities")
        # Utilities/u* threshold submenu
        self.menuUtilitiesUstar = QMenu(self.menuUtilities)
        self.menuUtilitiesUstar.setTitle("u* threshold")
        # Help menu
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setTitle("Help")
        # File menu items: menu actions for control files
        self.actionFileOpen = QAction(self)
        self.actionFileOpen.setText("Open")
        self.actionFileOpen.setShortcut('Ctrl+O')
        self.actionFileSave = QAction(self)
        self.actionFileSave.setText("Save")
        self.actionFileSave.setShortcut('Ctrl+S')
        self.actionFileSaveAs = QAction(self)
        self.actionFileSaveAs.setText("Save As...")
        self.actionFileSaveAs.setShortcut('Shift+Ctrl+S')
        # File/Convert submenu
        self.actionFileConvertnc2biomet = QAction(self)
        self.actionFileConvertnc2biomet.setText("nc to Biomet")
        self.actionFileConvertnc2xls = QAction(self)
        self.actionFileConvertnc2xls.setText("nc to Excel")
        self.actionFileConvertnc2reddyproc = QAction(self)
        self.actionFileConvertnc2reddyproc.setText("nc to REddyProc")
        # File menu item: split netCDF
        self.actionFileSplit = QAction(self)
        self.actionFileSplit.setText("Split")
        self.actionFileSplit.setShortcut('Ctrl+P')
        # File menu item: Quit
        self.actionFileQuit = QAction(self)
        self.actionFileQuit.setText("Quit")
        # the Vinod mod ...
        self.actionFileQuit.setShortcut('Ctrl+Q')
        # Edit menu items
        self.actionEditPreferences = QAction(self)
        self.actionEditPreferences.setText("Preferences...")
        # Run menu items
        self.actionRunCurrent = QAction(self)
        self.actionRunCurrent.setText("Current...")
        self.actionRunCurrent.setShortcut('Ctrl+R')
        self.actionStopCurrent = QAction(self)
        self.actionStopCurrent.setText("Stop run")
        # Plot menu items
        self.actionPlotFcVersusUstar = QAction(self)
        self.actionPlotFcVersusUstar.setText("Fco2 vs u*")
        self.actionPlotFingerprints = QAction(self)
        self.actionPlotFingerprints.setText("Fingerprints")
        self.actionPlotQuickCheck = QAction(self)
        self.actionPlotQuickCheck.setText("Summary")
        self.actionPlotTimeSeries = QAction(self)
        self.actionPlotTimeSeries.setText("Time series")
        self.actionPlotClosePlots = QAction(self)
        self.actionPlotClosePlots.setText("Close plots")
        # Utilities menu
        self.actionUtilitiesClimatology = QAction(self)
        self.actionUtilitiesClimatology.setText("Climatology")
        self.actionUtilitiesUstarCPDBarr = QAction(self)
        self.actionUtilitiesUstarCPDBarr.setText("CPD (Barr)")
        self.actionUtilitiesUstarCPDMcHugh = QAction(self)
        self.actionUtilitiesUstarCPDMcHugh.setText("CPD (McHugh)")
        self.actionUtilitiesUstarCPDMcNew = QAction(self)
        self.actionUtilitiesUstarCPDMcNew.setText("CPD (McNew)")
        self.actionUtilitiesUstarMPT = QAction(self)
        self.actionUtilitiesUstarMPT.setText("MPT")
        self.actionUtilitiesCFCheck = QAction(self)
        self.actionUtilitiesCFCheck.setText("CF check")
        # add the actions to the menus
        # File/Convert submenu
        self.menuFileConvert.addAction(self.actionFileConvertnc2xls)
        self.menuFileConvert.addAction(self.actionFileConvertnc2biomet)
        self.menuFileConvert.addAction(self.actionFileConvertnc2reddyproc)
        # File menu
        self.menuFile.addAction(self.actionFileOpen)
        self.menuFile.addAction(self.actionFileSave)
        self.menuFile.addAction(self.actionFileSaveAs)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.menuFileConvert.menuAction())
        self.menuFile.addAction(self.actionFileSplit)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionFileQuit)
        # Edit menu
        self.menuEdit.addAction(self.actionEditPreferences)
        # Run menu
        self.menuRun.addAction(self.actionRunCurrent)
        # Plot menu
        self.menuPlot.addAction(self.actionPlotFcVersusUstar)
        self.menuPlot.addAction(self.actionPlotFingerprints)
        self.menuPlot.addAction(self.actionPlotQuickCheck)
        self.menuPlot.addAction(self.actionPlotTimeSeries)
        self.menuPlot.addSeparator()
        self.menuPlot.addAction(self.actionPlotClosePlots)
        # Utilities/u* threshold submenu
        self.menuUtilitiesUstar.addAction(self.actionUtilitiesUstarCPDBarr)
        self.menuUtilitiesUstar.addAction(self.actionUtilitiesUstarCPDMcHugh)
        self.menuUtilitiesUstar.addAction(self.actionUtilitiesUstarCPDMcNew)
        self.menuUtilitiesUstar.addAction(self.actionUtilitiesUstarMPT)
        # Utilities menu
        self.menuUtilities.addAction(self.actionUtilitiesClimatology)
        self.menuUtilities.addAction(self.menuUtilitiesUstar.menuAction())
        self.menuUtilities.addAction(self.actionUtilitiesCFCheck)
        # add individual menus to menu bar
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuRun.menuAction())
        self.menubar.addAction(self.menuPlot.menuAction())
        self.menubar.addAction(self.menuUtilities.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        # create a tab bar
        self.tabs = QTabWidget(self)
        self.tabs.tab_index_all = 0
        self.tabs.tab_index_current = 0
        self.tabs.tab_dict = {}
        self.tabs.cfg_dict = {}
        # make the tabs closeable
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.closeTab)
        # add the text editor to the first tab
        self.tabs.addTab(textBox, "Log")
        self.tabs.tab_index_all = self.tabs.tab_index_all + 1
        # hide the tab close icon for the console tab
        self.tabs.tabBar().setTabButton(0, QTabBar.RightSide, None)
        # connect the tab-in-focus signal to the appropriate slot
        self.tabs.currentChanged[int].connect(self.tabSelected)

        # use VBoxLayout to position widgets so they resize with main window
        layout = QVBoxLayout()
        # add widgets to the layout
        layout.addWidget(self.menubar)
        layout.addWidget(self.tabs)
        self.setLayout(layout)
        self.setGeometry(50,50,800, 600)
        self.setWindowTitle(pfp_version)

        self.threadpool = QThreadPool()

        # Connect signals to slots
        # File menu actions
        self.actionFileConvertnc2biomet.triggered.connect(lambda:pfp_top_level.do_file_convert_nc2biomet(None, mode="standard"))
        self.actionFileConvertnc2xls.triggered.connect(pfp_top_level.do_file_convert_nc2xls)
        self.actionFileConvertnc2reddyproc.triggered.connect(lambda:pfp_top_level.do_file_convert_nc2reddyproc(None, mode="standard"))
        self.actionFileOpen.triggered.connect(self.file_open)
        self.actionFileSave.triggered.connect(self.file_save)
        self.actionFileSaveAs.triggered.connect(self.file_save_as)
        self.actionFileSplit.triggered.connect(pfp_top_level.do_file_split)
        self.actionFileQuit.triggered.connect(QApplication.quit)
        # Edit menu actions
        self.actionEditPreferences.triggered.connect(self.edit_preferences)
        # Run menu actions
        self.actionRunCurrent.triggered.connect(self.run_current)
        # Plot menu actions
        self.actionPlotFcVersusUstar.triggered.connect(pfp_top_level.do_plot_fcvsustar)
        self.actionPlotFingerprints.triggered.connect(pfp_top_level.do_plot_fingerprints)
        self.actionPlotQuickCheck.triggered.connect(pfp_top_level.do_plot_quickcheck)
        self.actionPlotTimeSeries.triggered.connect(pfp_top_level.do_plot_timeseries)
        self.actionPlotClosePlots.triggered.connect(pfp_top_level.do_plot_closeplots)
        # Utilities menu actions
        self.actionUtilitiesClimatology.triggered.connect(self.utilities_climatology)
        self.actionUtilitiesUstarCPDBarr.triggered.connect(self.utilities_ustar_cpd_barr)
        self.actionUtilitiesUstarCPDMcHugh.triggered.connect(self.utilities_ustar_cpd_mchugh_standard)
        self.actionUtilitiesUstarCPDMcNew.triggered.connect(self.utilities_ustar_cpd_mcnew)
        self.actionUtilitiesUstarMPT.triggered.connect(self.utilities_ustar_mpt)
        self.actionUtilitiesCFCheck.triggered.connect(self.utilities_cfcheck)
        # add the L4 GUI
        self.l4_ui = pfp_gui.pfp_l4_ui(self)
        # add the L5 GUI
        self.solo_gui = pfp_gui.solo_gui(self)

    def file_open(self, file_uri=None):
        """
        Purpose:
         Get a file name, figure out if it is a control file or a netCDF
         file and try to open it.
        """
        # local logger
        logger = logging.getLogger(name="pfp_log")
        # get the control file path
        if not file_uri:
            file_uri = QFileDialog.getOpenFileName(caption="Choose a file ...")[0]
            # check to see if file open was cancelled
            if len(str(file_uri)) == 0:
                return
        self.file_uri = str(file_uri)
        # read the contents of the control file
        logger.info(" Opening " + self.file_uri)
        # check to see if it is a control file
        file_open_success = False
        if not file_open_success:
            try:
                self.file = ConfigObj(self.file_uri, indent_type="    ", list_values=False,
                                      write_empty_values=True)
                file_open_success = True
            except Exception:
                # trying to open a netCDF file will throw UnicodeDecodeError
                # opening a ConfigObj file with syntax erros throws ConfigObjError
                exc_type, _, _ = sys.exc_info()
                if exc_type.__name__ == "UnicodeDecodeError":
                    # most likely a netCDF file
                    pass
                elif exc_type.__name__ == "ConfigObjError":
                    # syntax error in control file
                    msg = "Syntax error in control file, see below for line number"
                    logger.error(msg)
                    error_message = traceback.format_exc()
                    logger.error(error_message)
                    return
                else:
                    # unrecognised file type
                    msg = "File must be either a control file or a netCDF file"
                    logger.error(msg)
                    error_message = traceback.format_exc()
                    logger.error(error_message)
                    return
        # check to see if it is a netCDF file
        if not file_open_success:
            try:
                self.file = netCDF4.Dataset(self.file_uri, "r")
                file_open_success = True
            except:
                # unrecognised file type
                msg = "File must be either a control file or a netCDF file"
                logger.error(msg)
                error_message = traceback.format_exc()
                logger.error(error_message)
                return
        # do the business depending on the file type
        if isinstance(self.file, ConfigObj):
            # we are opening a control file
            self.open_control_file()
        elif isinstance(self.file, netCDF4._netCDF4.Dataset):
            # close the netCDF file
            self.file.close()
            # we are opening a netCDF file
            self.open_netcdf_file()
        else:
            # unrecognised file type
            msg = "File must be either a control file or a netCDF file"
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
        return

    def open_control_file(self):
        logger = logging.getLogger(name="pfp_log")
        self.cfg = self.file
        # check to see if the processing level is defined in the control file
        if "level" not in self.cfg:
            # if not, then sniff the control file to see what it is
            self.cfg["level"] = self.get_cf_level()
            # and save the control file
            self.cfg.write()
        # create a QtTreeView to edit the control file
        if self.cfg["level"] in ["L1"]:
            # update control file to new syntax
            if not pfp_compliance.l1_update_controlfile(self.cfg): return
            # put the GUI for editing the L1 control file in a new tab
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L1(self)
            # get the control file data from the L1 edit GUI
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["L2"]:
            if not pfp_compliance.l2_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L2(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["L3"]:
            if not pfp_compliance.l3_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L3(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["concatenate"]:
            if not pfp_compliance.concatenate_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_concatenate(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["climatology"]:
            if not pfp_compliance.climatology_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_climatology(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["cpd_barr"]:
            if not pfp_compliance.cpd_barr_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_cpd_barr(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["cpd_mchugh"]:
            if not pfp_compliance.cpd_mchugh_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_cpd_mchugh(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["cpd_mcnew"]:
            if not pfp_compliance.cpd_mcnew_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_cpd_mcnew(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["mpt"]:
            if not pfp_compliance.mpt_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_mpt(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["L4"]:
            if not pfp_compliance.l4_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L4(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["L5"]:
            if not pfp_compliance.l5_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L5(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["L6"]:
            if not pfp_compliance.l6_update_controlfile(self.cfg): return
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_L6(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["nc2csv_biomet"]:
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_nc2csv_biomet(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["nc2csv_ecostress"]:
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_nc2csv_ecostress(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["nc2csv_fluxnet"]:
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_nc2csv_fluxnet(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["nc2csv_reddyproc"]:
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_nc2csv_reddyproc(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        elif self.cfg["level"] in ["batch"]:
            self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.edit_cfg_batch(self)
            self.tabs.cfg_dict[self.tabs.tab_index_all] = self.tabs.tab_dict[self.tabs.tab_index_all].get_data_from_model()
        else:
            logger.error(" Unrecognised control file type: " + self.cfg["level"])
            return
        # file type and uri
        self.tabs.cfg_dict[self.tabs.tab_index_all]["file"] = self.file
        self.tabs.cfg_dict[self.tabs.tab_index_all]["file_uri"] = self.file_uri
        # add a tab for the control file
        self.tabs.addTab(self.tabs.tab_dict[self.tabs.tab_index_all], os.path.basename(str(self.file_uri)))
        self.tabs.setCurrentIndex(self.tabs.tab_index_all)
        if self.tabs.tab_dict[self.tabs.tab_index_all].cfg_changed:
            self.update_tab_text()
        self.tabs.tab_index_all = self.tabs.tab_index_all + 1
        return

    def open_netcdf_file(self):
        # read the netCDF file to a data structure
        self.ds = pfp_io.nc_read_series(self.file_uri)
        if self.ds.returncodes["value"] != 0:
            return
        file_path_parts = os.path.split(self.file_uri)
        self.file_path = file_path_parts[0]
        self.file_name = file_path_parts[1]
        # display the netcdf file in the GUI
        self.tabs.tab_dict[self.tabs.tab_index_all] = pfp_gui.file_explore(self)
        # return if something went wrong
        if self.tabs.tab_dict[self.tabs.tab_index_all].ds.returncodes["value"] != 0:
            return
        # create the cfg dictionary
        self.tabs.cfg_dict[self.tabs.tab_index_all] = {}
        # put the file path into the tabs dictionary
        self.tabs.cfg_dict[self.tabs.tab_index_all]["file"] = self.file
        self.tabs.cfg_dict[self.tabs.tab_index_all]["file_uri"] = self.file_uri
        # add a tab for the netCDF file contents
        self.tabs.addTab(self.tabs.tab_dict[self.tabs.tab_index_all], self.file_name)
        self.tabs.setCurrentIndex(self.tabs.tab_index_all)
        self.tabs.tab_index_all = self.tabs.tab_index_all + 1
        return

    def get_cf_level(self):
        """ Sniff the control file to find out it's type."""
        logger = logging.getLogger(name="pfp_log")
        self.cfg["level"] = ""
        # check for L1
        if self.check_cfg_L1():
            logger.info(" L1 control file detected")
            self.cfg["level"] = "L1"
        # check for L2
        elif self.check_cfg_L2():
            logger.info(" L2 control file detected")
            self.cfg["level"] = "L2"
        # check for L3
        elif self.check_cfg_L3():
            logger.info(" L3 control file detected")
            self.cfg["level"] = "L3"
        # check for concatenate
        elif self.check_cfg_concatenate():
            logger.info(" Concatenate control file detected")
            self.cfg["level"] = "concatenate"
        # check for L4
        elif self.check_cfg_L4():
            logger.info(" L4 control file detected")
            self.cfg["level"] = "L4"
        # check for L5
        elif self.check_cfg_L5():
            logger.info(" L5 control file detected")
            self.cfg["level"] = "L5"
        # check for L6
        elif self.check_cfg_L6():
            logger.info(" L6 control file detected")
            self.cfg["level"] = "L6"
        else:
            logger.info(" Unable to detect level, enter manually ...")
            text, ok = QInputDialog.getText(self, 'Processing level', 'Enter the processing level:')
            if ok:
                self.cfg["level"] = text
        return self.cfg["level"]

    def check_cfg_L1(self):
        """ Return true if a control file is an L1 file."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            # remove the common sections
            common_sections = ["level", "controlfile_name", "Files", "Global", "Output",
                               "Plots", "General", "Options", "Soil", "Massman", "GUI"]
            for section in list(self.cfg.keys()):
                if section in common_sections:
                    cfg_sections.remove(section)
            # now loop over the remaining sections
            for section in cfg_sections:
                subsections = list(self.cfg[section].keys())
                for subsection in subsections:
                    if "Attr" in list(self.cfg[section][subsection].keys()):
                        result = True
                        break
        except:
            result = False
        return result

    def check_cfg_L2(self):
        """ Return true if a control file is an L2 file."""
        result = False
        try:
            got_sections = False
            cfg_sections = list(self.cfg.keys())
            if (("Files" in cfg_sections) and
                ("Variables" in cfg_sections)):
                got_sections = True
            # loop over [Variables] sections
            got_qc = False
            qc_list = ["RangeCheck", "DiurnalCheck", "ExcludeDates", "DependencyCheck", "UpperCheck",
                       "LowerCheck", "ExcludeHours", "Linear", "CorrectWindDirection"]
            for section in ["Variables"]:
                subsections = list(self.cfg[section].keys())
                for subsection in subsections:
                    for qc in qc_list:
                        if qc in list(self.cfg[section][subsection].keys()):
                            got_qc = True
                            break
            # final check
            if got_sections and got_qc and not self.check_cfg_L3() and not self.check_cfg_L4():
                result = True
        except:
            result = False
        return result

    def check_cfg_L3(self):
        """ Return true if a control file is an L3 file."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            if ((("General" in cfg_sections) or
                ("Soil" in cfg_sections) or
                ("Massman" in cfg_sections)) and
                ("Options" in cfg_sections)):
                result = True
        except:
            result = False
        return result

    def check_cfg_concatenate(self):
        """ Return true if control file is concatenation."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            if "Files" in cfg_sections:
                if (("Out" in list(self.cfg["Files"].keys())) and
                    ("In" in list(self.cfg["Files"].keys()))):
                    result = True
        except:
            result = False
        return result

    def check_cfg_L4(self):
        """ Return true if control file is L4."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            for section in cfg_sections:
                if section in ["Variables", "Drivers", "Fluxes"]:
                    subsections = list(self.cfg[section].keys())
                    for subsection in subsections:
                        if (("GapFillFromAlternate" in list(self.cfg[section][subsection].keys())) or
                            ("GapFillFromClimatology" in list(self.cfg[section][subsection].keys()))):
                            result = True
                            break
        except:
            result = False
        return result

    def check_cfg_L5(self):
        """ Return true if control file is L5."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            for section in cfg_sections:
                if section in ["Variables", "Drivers", "Fluxes"]:
                    subsections = list(self.cfg[section].keys())
                    for subsection in subsections:
                        if (("GapFillUsingSOLO" in list(self.cfg[section][subsection].keys())) or
                            ("GapFillUsingMDS" in list(self.cfg[section][subsection].keys()))):
                            result = True
                            break
        except:
            result = False
        return result

    def check_cfg_L6(self):
        """ Return true if control file is L6."""
        result = False
        try:
            cfg_sections = list(self.cfg.keys())
            if ("EcosystemRespiration" in cfg_sections or
                "NetEcosystemExchange" in cfg_sections or
                "GrossPrimaryProductivity" in cfg_sections):
                result = True
        except:
            result = False
        return result

    def direct_run(self):
        """ Placeholder until full implementation done."""
        logger = logging.getLogger(name="pfp_log")
        msg = " Open control file and use 'Run/Current ...'"
        logger.warning(msg)
        return

    def file_save(self):
        """ Save the current file."""
        logger = logging.getLogger(name="pfp_log")
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        if isinstance(self.tabs.cfg_dict[tab_index_current]["file"],
                      ConfigObj):
            # we are opening a control file
            self.save_control_file()
        elif isinstance(self.tabs.cfg_dict[tab_index_current]["file"],
                        netCDF4._netCDF4.Dataset):
            # we are opening a netCDF file
            self.save_netcdf_file()
        else:
            # unrecognised file type
            msg = "File must be either a control file or a netCDF file"
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
        return

    def save_control_file(self):
        """ Save the current tab as a control file."""
        logger = logging.getLogger(name="pfp_log")
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        # get the control file name
        cfg_filename = self.tabs.cfg_dict[tab_index_current]["file_uri"]
        # get the updated control file data

        cfg = self.tabs.tab_dict[tab_index_current].get_data_from_model()

        # check to make sure we are not overwriting the template version
        if "template" not in cfg_filename:
            # set the control file name
            cfg.filename = cfg_filename
        else:
            msg = " You are trying to write to the template folder.\n"
            msg = msg + "Please save this control file to a different location."
            msgbox = pfp_gui.myMessageBox(msg)
            # put up a "Save as ..." dialog
            cfg_filename = QFileDialog.getSaveFileName(self, "Save as ...")[0]
            # return without doing anything if cancel used
            if len(str(cfg_filename)) == 0:
                return
            # set the control file name
            cfg.filename = str(cfg_filename)
        # write the control file
        msg = " Saving " + cfg.filename
        logger.info(msg)
        cfg.write()
        # remove the asterisk in the tab text
        tab_text = str(self.tabs.tabText(tab_index_current))
        self.tabs.setTabText(self.tabs.tab_index_current, tab_text.replace("*",""))
        # reset the cfg changed logical to false
        self.tabs.tab_dict[tab_index_current].cfg_changed = False
        return

    def save_netcdf_file(self):
        """Save the current tab as a netCDF file."""
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        # get the netCDF file name
        nc_file_uri = self.tabs.cfg_dict[tab_index_current]["file_uri"]
        # get the updated control file data
        ds = self.tabs.tab_dict[tab_index_current].get_data_from_model()
        # write the data structure to file
        pfp_io.NetCDFWrite(nc_file_uri, ds)
        # remove the asterisk in the tab text
        tab_text = str(self.tabs.tabText(tab_index_current))
        self.tabs.setTabText(self.tabs.tab_index_current, tab_text.replace("*",""))
        # reset the cfg changed logical to false
        self.tabs.tab_dict[tab_index_current].cfg_changed = False
        return

    def file_save_as(self):
        """Save a file with a different name."""
        logger = logging.getLogger(name="pfp_log")
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        # put up a "Save as ..." dialog
        file_uri = self.tabs.cfg_dict[tab_index_current]["file_uri"]
        file_uri = QFileDialog.getSaveFileName(self, "Save as ...", file_uri)[0]
        # return without doing anything if cancel used
        if len(str(file_uri)) == 0:
            return
        self.tabs.cfg_dict[tab_index_current]["file_uri"] = file_uri
        if isinstance(self.tabs.cfg_dict[tab_index_current]["file"],
                      ConfigObj):
            # we are opening a control file
            self.save_as_control_file()
        elif isinstance(self.tabs.cfg_dict[tab_index_current]["file"],
                        netCDF4._netCDF4.Dataset):
            # we are opening a netCDF file
            self.save_as_netcdf_file()
        else:
            # unrecognised file type
            msg = "File must be either a control file or a netCDF file"
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
        return

    def save_as_control_file(self):
        """ Save the current tab with a different name."""
        logger = logging.getLogger(name="pfp_log")
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        file_uri = self.tabs.cfg_dict[tab_index_current]["file_uri"]
        # get the updated control file data
        cfg = self.tabs.tab_dict[tab_index_current].get_data_from_model()
        # set the control file name
        cfg.filename = str(file_uri)
        # write the control file
        logger.info(" Saving " + cfg.filename)
        cfg.write()
        # update the control file name
        self.tabs.cfg_dict[tab_index_current]["file_uri"] = cfg.filename
        # update the tab text
        self.tabs.setTabText(tab_index_current, os.path.basename(str(file_uri)))
        # reset the cfg changed logical to false
        self.tabs.tab_dict[tab_index_current].cfg_changed = False
        return

    def save_as_netcdf_file(self):
        """ Save the current tab with a different name."""
        # get the current tab index
        tab_index_current = self.tabs.tab_index_current
        # get the netCDF file name
        file_uri = self.tabs.cfg_dict[tab_index_current]["file_uri"]
        # get the updated control file data
        ds = self.tabs.tab_dict[tab_index_current].get_data_from_model()
        # write the data structure to file
        pfp_io.NetCDFWrite(file_uri, ds)
        # update the tab text
        self.tabs.setTabText(tab_index_current, os.path.basename(str(file_uri)))
        # reset the cfg changed logical to false
        self.tabs.tab_dict[tab_index_current].cfg_changed = False
        return

    def edit_preferences(self):
        print("Edit/Preferences goes here")
        pass

    def tabSelected(self, arg=None):
        self.tabs.tab_index_current = arg

    def run_current(self):
        self.stop_flag = False
        # save the current tab index
        logger = logging.getLogger(name="pfp_log")
        tab_index_current = self.tabs.tab_index_current
        if tab_index_current == 0:
            msg = " No control file selected ..."
            logger.warning(msg)
            return
        # disable the Run/Current menu option
        self.actionRunCurrent.setDisabled(True)
        # get the updated control file data
        self.cfg = self.tabs.tab_dict[tab_index_current].get_data_from_model()
        # set the focus back to the log tab
        self.tabs.setCurrentIndex(0)
        # call the appropriate processing routine depending on the level
        self.tabs.tab_index_running = tab_index_current
        if self.tabs.cfg_dict[tab_index_current]["level"] == "batch":
            # check the L1 control file to see if it is OK to run
            if not pfp_compliance.check_batch_controlfile(self.cfg): return
            # add stop to run menu
            self.menuRun.addAction(self.actionStopCurrent)
            self.actionStopCurrent.triggered.connect(self.stop_current)
            # get a worker thread
            worker = pfp_threading.Worker(pfp_top_level.do_run_batch, self)
            # start the worker
            self.threadpool.start(worker)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L1":
            # check the L1 control file to see if it is OK to run
            if pfp_compliance.check_l1_controlfile(self.cfg):
                pfp_top_level.do_run_l1(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L2":
            pfp_top_level.do_run_l2(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L3":
            pfp_top_level.do_run_l3(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "concatenate":
            pfp_top_level.do_file_concatenate(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "climatology":
            self.utilities_climatology(cfg=self.cfg, mode="custom")
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "cpd_barr":
            self.utilities_ustar_cpd_barr(cfg=self.cfg, mode="custom")
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "cpd_mchugh":
            self.utilities_ustar_cpd_mchugh_custom()
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "cpd_mcnew":
            self.utilities_ustar_cpd_mcnew(cfg=self.cfg, mode="custom")
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "mpt":
            pfp_top_level.do_utilities_ustar_mpt(cfg=self.cfg, mode="custom")
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L4":
            pfp_top_level.do_run_l4(self)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L5":
            pfp_top_level.do_run_l5(self)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "L6":
            pfp_top_level.do_run_l6(self)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "nc2csv_biomet":
            pfp_top_level.do_file_convert_nc2biomet(self.cfg, mode="custom")
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "nc2csv_ecostress":
            pfp_top_level.do_file_convert_nc2ecostress(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "nc2csv_fluxnet":
            pfp_top_level.do_file_convert_nc2fluxnet(self.cfg)
        elif self.tabs.cfg_dict[tab_index_current]["level"] == "nc2csv_reddyproc":
            pfp_top_level.do_file_convert_nc2reddyproc(self.cfg, mode="custom")
        else:
            logger.error("Level not implemented yet ...")
        # enable the Run/Current menu option
        self.actionRunCurrent.setDisabled(False)
        return

    def stop_current(self):
        # put up a message box, continue or quit
        msg = "Do you want to quit processing?"
        result = pfp_gui.MsgBox_ContinueOrQuit(msg, title="Quit processing?")
        if result.clickedButton().text() == "Quit":
            self.stop_flag = True
            msg = "Processing will stop when this level is finished"
            result = pfp_gui.myMessageBox(msg)
        return

    def closeTab (self, currentIndex):
        """ Close the selected tab."""
        # check to see if the tab contents have been saved
        tab_text = str(self.tabs.tabText(currentIndex))
        if "*" in tab_text:
            msg = "Save control file?"
            reply = QMessageBox.question(self, 'Message', msg,
                                               QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.file_save()
        # get the current tab from its index
        currentQWidget = self.tabs.widget(currentIndex)
        # delete the tab
        currentQWidget.deleteLater()
        self.tabs.removeTab(currentIndex)
        # remove the corresponding entry in cfg_dict
        self.tabs.cfg_dict.pop(currentIndex)
        # and renumber the keys
        for n in list(self.tabs.cfg_dict.keys()):
            if n > currentIndex:
                self.tabs.cfg_dict[n-1] = self.tabs.cfg_dict.pop(n)
        # remove the corresponding entry in tab_dict
        self.tabs.tab_dict.pop(currentIndex)
        # and renumber the keys
        for n in list(self.tabs.tab_dict.keys()):
            if n > currentIndex:
                self.tabs.tab_dict[n-1] = self.tabs.tab_dict.pop(n)
        # decrement the tab index
        self.tabs.tab_index_all = self.tabs.tab_index_all - 1
        return

    def update_tab_text(self):
        """ Add an asterisk to the tab title text to indicate tab contents have changed."""
        # add an asterisk to the tab text to indicate the tab contents have changed
        tab_text = str(self.tabs.tabText(self.tabs.tab_index_current))
        if "*" not in tab_text:
            self.tabs.setTabText(self.tabs.tab_index_current, tab_text+"*")

    def utilities_climatology(self, cfg=None, mode="standard"):
        logger = logging.getLogger(name="pfp_log")
        if mode == "standard":
            nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
            if not os.path.exists(nc_file_uri):
                logger.info( " Climatology: no netCDF file chosen")
                return
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_climatology,
                                          nc_file_uri=nc_file_uri, mode=mode)
        else:
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_climatology,
                                          cfg=cfg, mode=mode)
        self.threadpool.start(worker)

    def utilities_ustar_cpd_barr(self, cfg=None, mode="standard"):
        logger = logging.getLogger(name="pfp_log")
        if mode == "standard":
            nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
            if not os.path.exists(nc_file_uri):
                logger.info( " CPD (Barr): no netCDF file chosen")
                return
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_barr,
                                          nc_file_uri=nc_file_uri, mode=mode)
        else:
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_barr,
                                          cfg=cfg, mode=mode)
        self.threadpool.start(worker)

    def utilities_ustar_cpd_mchugh_standard(self):
        # disable the Run/Current menu option
        self.actionRunCurrent.setDisabled(True)
        nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
        if not os.path.exists(nc_file_uri):
            logger.info( " CPD (McHugh): no netCDF file chosen")
            return
        worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_mchugh,
                                      self, nc_file_uri=nc_file_uri, mode="standard")
        self.threadpool.start(worker)
        return
    def utilities_ustar_cpd_mchugh_custom(self):
        # disable the Run/Current menu option
        self.actionRunCurrent.setDisabled(True)
        worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_mchugh,
                                      self, cfg=self.cfg, mode="custom")
        self.threadpool.start(worker)
        return
    def utilities_ustar_cpd_mcnew(self, cfg=None, mode="standard"):
        logger = logging.getLogger(name="pfp_log")
        if mode == "standard":
            nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
            if not os.path.exists(nc_file_uri):
                logger.info( " CPD (McNew): no netCDF file chosen")
                return
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_mcnew,
                                          nc_file_uri=nc_file_uri, mode=mode)
        else:
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_cpd_mcnew,
                                          cfg=cfg, mode=mode)
        self.threadpool.start(worker)

    def utilities_ustar_mpt(self, cfg=None, mode="standard"):
        if mode == "standard":
            nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
            if not os.path.exists(nc_file_uri):
                logger.info( " CPD (MPT): no netCDF file chosen")
                return
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_mpt,
                                          nc_file_uri=nc_file_uri, mode=mode)
        else:
            worker = pfp_threading.Worker(pfp_top_level.do_utilities_ustar_mpt,
                                          cfg=cfg, mode=mode)
        self.threadpool.start(worker)

    def utilities_cfcheck(self):
        nc_file_uri = pfp_io.get_filename_dialog(title="Choose a netCDF file", ext="*.nc")
        if not os.path.exists(nc_file_uri):
            logger.info( " CF check: no netCDF file chosen")
            return
        worker = pfp_threading.Worker(pfp_top_level.do_utilities_cfcheck,
                                      nc_file_uri)
        self.threadpool.start(worker)

if (__name__ == '__main__'):
    # get the application name and version
    pfp_version = cfg.version_name + " " + cfg.version_number
    app = QApplication(["PyFluxPro"])
    # get a text editor for the log window
    textBox = QTextEdit()
    textBox.setReadOnly(True)
    logger.consoleHandler.sigLog.connect(textBox.append)
    # instance the main GUI
    ui = pfp_main_ui(pfp_version, textBox)
    ui.show()
    #pfp_compliance.check_executables()
    app.exec_()
    del ui
