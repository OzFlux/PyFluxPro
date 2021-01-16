version_name = "PyFluxPro"
version_number = "V3.0.6"
# V3.0.6 - 13th January 2021
#        - patch code for getting the height of the CO2 sensor
#          for calculating single point storage
# V3.0.5 - 25th November 2020
#        - implement plotting at L2 and L3 in batch mode
# V3.0.4 - 20th November 2020
#        - PyInstaller version
#        - corrected syntax for importing PFP modules
#        - added pfp_utils.get_base_path() to return correct appilcation
#          path whether running native Python version or PyInstaller version
# V3.0.3 - 25th October 2020 port of changes in PFP2.7 "another_fine_mess" branch
#        - major changes (to be documented later)
#        - 6th August 2020 (long overdue version change)
#        - check GitHub for what has been done since 5th August 2019
#        - most recent changes are:
#          - implementation of pfp_gui.TreeView() subclassed from QTreeView
#            to get the drag and drop behaviour we want
#          - rename Eva's second ustar filter method to "ustar (EvGb)"
# V3.0.2 - 14th Sep 2020
#          introduced changes in PFP2.7 by Peter made in Aug/Sep
#          to python3 before he takes over
#          e.g. new feature "explore"
# V3.0.1 - 21st May 2020
#          converted PFP2.7 to python3
#          using 2to3 and manual (esp. integer division)
#          conversion
# V1.0.1 - 5th August 2019
#        - remove redundant SOLO code at L6
#        - implement L1 to L4 batch processing
# V1.0.0 - 16th June 2019
#        - MAJOR CHANGES
#          - implement detection and filling of long gaps
#            at L5
#          - implement consistent way of handling program
#            settings and options at L4, L5 and L6
#          - move from PyQt4 to PyQt5
#          - many other changes along the way
# V0.9.2 - 7th October 2018
#        - implemented output to ECOSTRESS CSV file
#        - implemented L5 and L6 processing
#        - cleaned up the GUI routines in pfp_gui.py
# V0.9.1 - updated to PyFluxPro V0.1.5
#        - working steps at this stage are:
#          - L1 to L4
#          - concatenation
#          - climatology
# V0.9.0 - fork from main PyFluxPro at V0.1.1
#        - implement integrated GUI
# V0.1.1 - cumulative changes up to 13/07/2017
#        - major change to pfp_utils.CreateSeries()
#          - deprecated "FList=" argument
#          - made all arguments compulsory
#        - QC flag generated for each series immediately
#          prior to pfp_utils.CreateSeries() based solely
#          on series data mask
# V0.1.0 - copy of OzFluxQC V2.9.6f
