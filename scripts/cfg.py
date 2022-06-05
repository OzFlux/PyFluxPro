version_name = "PyFluxPro"
version_number = "V3.4.2"
# V3.4.2 - June 2022
#        - added Cacilia's windrose plotting code
#          - standard (via menu)
#          - custom (via control file)
# V3.4.1 - May 2022
#        - fixed missing data after 2D interpolation in pfp_clim.climatology()
# V3.4.0 - April 2022
#        - replaced Ian McHugh's original partitioning code with his "new"
#          (2019) partitioning code
#        - tested with Calperum, Cumberland Plain and Loxton
#          - little difference for Lloyd-Taylor
#          - Lasslop et al now much closer to Lloyd-Taylor
#        - changes to editing of L6 control file
# V3.3.5 - March 2022
#        - added Linear to L1, works the same as Linear at L2
# V3.3.4 - March 2022
#        - rewrite pfp_cpd_mchugh.fit() to use ndarray rather than pandas df
# V3.3.3 - February 2022
#        - fix standard_name for Fg_Av in pfp_ts.CorrectFgForStorage()
#        - add statistic_type to Fco2_single in pfp_ts.CalculateFco2StorageSinglePoint()
#        - fixed data_link and fluxnet_id in utilities/cleanup_netcdf_files.py
# V3.3.2 - January 2022
#        - implement EasyFlux-DL variables in check_l1_controlfile.txt
# V3.3.1 - November 2021, post workshop bug fixes
#        - mainy cleaning up L1 checks before running
#        - trap frac ==> m^3/m^3 for RH
#        - make sure statistic_type added for variances
# V3.3.0 - mid-June 2021
#        - changes required for DSA compliance phase 1
#          - deprecated "group_name" and "serial_number" variable attributes
#          - changed code to deal with no "units" attribute
#          - "valid_range" variable attribute only set by QC routines
#          - sundry other changes
# V3.2.0 - 8th June 2021
#        - multithreading for batch processing
#          - rewrote initialisation of logger as a class (Whooot, go you OO guru!)
#          - rewrote logging to application "Log" window
#          - created new module pfp_threading
#        - multithreading for all items under "Utilities"
#        - multithreading for all utilities run using "Run/Current"
#        - added the ability to change netCDF variable names
# V3.1.7 - 3rd June 2021
#        - reading of AmeriFlux & FluxNet CSV files at L1
#        - fixed bug when height of CO2 sensor could not be determined
#        - set interpolation to "none" in pfp_plot.plot_explore_fingerprints()
# V3.1.6 - 29th May 2021
#        - implement Ian's "new" CPD code
# V3.1.5 - 21st May 2021
#        - File/Open, File/Save, File/Save as now handle control and netCDF files
#        - updated QC flags at L5 and L6 to use new convention (first digit = level)
# V3.1.4 - 20th May 2021
#        - updated utilities/cleanup_netcdf_files.py
#        - minor tweaks to get batch mode working
# V3.1.3 - 11th May 2021
#        - rewrite of pfp_ts.InterpolateOverMissing()
# V3.1.2 - 27th April 2021
#        - fixed bug in pfp_ck.do_uppercheck()
#        - trap cancel when choosing file under File/Explore
#        - added plotting of fingerprint to File/Explore right-click contect menu
# V3.1.1 - 19th February 2021
#        - added reading of CSV files
#        - basic editing L1 control files
# V3.1.0 - 16th February 2021
#        - replace xlrd with pandas.read_excel() after reading
#          of XLSX files was removed from xlrd at V1.2.0
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
