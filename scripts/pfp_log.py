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
def init_logger(logger_name, log_file_name, to_file=True, to_screen=False):
    """
    Purpose:
     Returns a logger object.
    Usage:
     logger = pfp_log.init_logger()
    Author: PRI with acknowledgement to James Cleverly
    Date: September 2016
    """
    # create formatter
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
    # create the logger and set the level
    logger = logging.getLogger(name=logger_name)
    logger.setLevel(logging.DEBUG)
    if to_file:
        # create file handler for all messages
        fh1 = logging.FileHandler(log_file_name)
        fh1.setLevel(logging.DEBUG)
        fh1.setFormatter(formatter)
        # add the file handler to the logger
        logger.addHandler(fh1)
        # set up a separate file for errors
        ext = os.path.splitext(log_file_name)[1]
        error_file_name = log_file_name.replace(ext, ".errors")
        fh2 = logging.FileHandler(error_file_name)
        fh2.setLevel(logging.ERROR)
        fh2.setFormatter(formatter)
        logger.addHandler(fh2)
        # ... and a separate file for warnings
        ext = os.path.splitext(log_file_name)[1]
        warning_file_name = log_file_name.replace(ext, ".warnings")
        fh3 = logging.FileHandler(warning_file_name)
        fh3.setLevel(logging.WARNING)
        fh3.setFormatter(formatter)
        logger.addHandler(fh3)
    if to_screen:
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        console.setLevel(logging.DEBUG)
        logger.addHandler(console)
    return logger
def change_logger_filename(logger_name, new_file_name):
    # get the logger
    logger = logging.getLogger(name=logger_name)
    # remove the existing file handlers
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            old_log_path = handler.baseFilename
            old_log_level = handler.level
            old_log_formatter = handler.formatter
            logger.removeHandler(handler)
            old_dir_name = os.path.dirname(os.path.abspath(old_log_path))
            old_base_name = os.path.basename(os.path.abspath(old_log_path))
            old_file_name = os.path.splitext(old_base_name)[0]
            new_base_name = old_base_name.replace(old_file_name, new_file_name)
            new_log_path = os.path.join(old_dir_name, new_base_name)
            fh = logging.FileHandler(new_log_path)
            fh.setLevel(old_log_level)
            fh.setFormatter(old_log_formatter)
            logger.addHandler(fh)
    return logger
def debug_function_enter(function_name):
    msg = "   Entering " + function_name
    logger.debug(msg)
    return
def debug_function_leave(function_name):
    msg = "   Leaving " + function_name
    logger.debug(msg)
    return
def disable_console_log(logger_name):
    console_handler = None
    logger = logging.getLogger(name=logger_name)
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            console_handler = handler
            logger.removeHandler(handler)
    return console_handler
def enable_console_log(logger_name, console_handler):
    if console_handler is None:
        return
    logger = logging.getLogger(name=logger_name)
    logger.addHandler(console_handler)
    return
def get_batch_log_path(log_path):
    if not os.path.isdir(log_path):
        os.mkdir(log_path)
    batch_log_path = os.path.join(log_path, "batch")
    if not os.path.isdir(batch_log_path):
        os.mkdir(batch_log_path)
    now = datetime.datetime.now()
    batch_log_now_path = os.path.join(batch_log_path, now.strftime("%Y%m%d%H%M"))
    if not os.path.isdir(batch_log_now_path):
        os.mkdir(batch_log_now_path)
    return batch_log_now_path
