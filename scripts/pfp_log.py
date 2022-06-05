# standard modules
import datetime
import logging
import os
# 3rd party modules
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

logger = logging.getLogger("pfp_log")

class ConsoleWindowLogHandler(logging.Handler, QObject):
    sigLog = pyqtSignal(str)
    def __init__(self):
        logging.Handler.__init__(self)
        QObject.__init__(self)

    def emit(self, logRecord):
        message = self.format(logRecord)
        self.sigLog.emit(message)
        # not needed when using threading but retained for those
        # proceses where threading is not used so that log window
        # updates whenever a message is emitted
        QApplication.processEvents()

class CreateLogger(logging.getLoggerClass()):
    """
    Purpose:
     Returns a logger object.
    Usage:
     logger = CreateLogger()
    Author: PRI with acknowledgement to James Cleverly
    Date: September 2016
    Mods: June 2021 - rewritten to be a class
    """
    def __init__(self, logger_name, to_log_window=True, to_screen=False, to_file=True, log_file_name="pfp_log"):
        super(CreateLogger, self).__init__(logger_name)
        # create formatter
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
        logger_level = logging.INFO
        # create the logger and set the level
        self.logger = logging.getLogger(name=logger_name)
        # prevent this logger propagating messages to the root logger
        # NOTE: siphon use logging.basicConfig() which sets up a root logger
        #       without the lne below, PFP logging messages appear in the
        #       terminal window
        self.logger.propagate = False
        self.logger.setLevel(logger_level)
        # default is to log to the application "Log" window
        if to_log_window:
            self.consoleHandler = ConsoleWindowLogHandler()
            self.consoleHandler.setFormatter(formatter)
            self.logger.addHandler(self.consoleHandler)
        if to_file:
            # create file handler for all messages
            fh1 = logging.FileHandler(log_file_name)
            fh1.setLevel(logger_level)
            fh1.setFormatter(formatter)
            # add the file handler to the logger
            self.logger.addHandler(fh1)
            # set up a separate file for errors
            ext = os.path.splitext(log_file_name)[1]
            error_file_name = log_file_name.replace(ext, ".errors")
            fh2 = logging.FileHandler(error_file_name)
            fh2.setLevel(logging.ERROR)
            fh2.setFormatter(formatter)
            self.logger.addHandler(fh2)
            # ... and a separate file for warnings
            ext = os.path.splitext(log_file_name)[1]
            warning_file_name = log_file_name.replace(ext, ".warnings")
            fh3 = logging.FileHandler(warning_file_name)
            fh3.setLevel(logging.WARNING)
            fh3.setFormatter(formatter)
            self.logger.addHandler(fh3)
        if to_screen:
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logger_level)
            self.logger.addHandler(console)
def debug_function_enter(function_name):
    msg = "   Entering " + function_name
    logger.debug(msg)
    return
def debug_function_leave(function_name):
    msg = "   Leaving " + function_name
    logger.debug(msg)
    return
