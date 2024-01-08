import os
import sys
from configobj import ConfigObj
cwd = os.getcwd()
sys.path.insert(0, cwd[:cwd.index("utilities")-1])
import scripts.pfp_io as pfp_io
import scripts.pfp_rp as pfp_rp

cfg_uri = "/home/pisaac/PyFluxPro/controlfiles/OzFlux/Tower/L6/AliceSpringsMulga.txt"
cfg = ConfigObj(cfg_uri, indent_type="    ", list_values=False, write_empty_values=True)
if "Options" not in cfg:
    cfg["Options"] = {}
cfg["Options"]["call_mode"] = "interactive"
cfg["Options"]["show_plots"] = "Yes"

nc_uri = os.path.join(cfg["Files"]["file_path"], cfg["Files"]["out_filename"])
ds6 = pfp_io.NetCDFRead(nc_uri)

l6_info = pfp_rp.ParseL6ControlFile(cfg, ds6)

pfp_rp.L6_summary(ds6, l6_info)
