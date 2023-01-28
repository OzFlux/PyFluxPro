import copy

class DataStructure(object):
    """
    This is the main data structure used throughout PyFluxPro.

    A data structure is used to hold all of the data and metadata used by
    PyFluxPro during data procssing.  It is basically a representation in
    Python of a netCDF file.

    Creating a data structure:
     A data structure can be created directly using;
      ds = pfp_classes.Datastructure()

    Attributes of a data structure:
     ds.root["Attributes"] - a dictionary containing the global attributes read from
                           or written to a netCDF file.
     ds.root["Variables"] - a dictionary containing the variables read from or written to a
                            netCDF file.

    Variables in a data structure:
     Each variable in the ds.root["Variables"] dictionary is another dictionary that
     contains the data, the quality control flag and the attributes.
     For example, the air temperature variable (Ta) is stored in the data
     structure as follows;
      ds.root["Variables"]["Ta"]["Data"] - a 1D numpy array containing the data values
      ds.root["Variables"]["Ta"]["Flag"] - a 1D numpy array containing the quality
                                control flag values
      ds.root["Variables"]["Ta"]["Attr"] - a dictionary containing the variable attributes

    Reading a data structure from file:
     A data structure is returned when a netCDF file is read e.g.;
      ds = pfp_io.NetCDFRead(file_name)
      where file_name is the name of the netCDF file.

    Writing a data structure to file:
     A data structure can be written to a netCDF file using;
      pfp_io.NetCDFWrite(file_name, ds)
      where file_name is the name of the netCDF file to be created
            ds is the data structure to be written to the file

    Accessing variables in a data structure:
     Variables in a data structure can be accessed using the GetVariable
     helper function as follows;
      Ta = pfp_utils.GetVariable(ds, "Ta")
     where "Ta" is the label of the variable.
           ds is the data structure containing the variable.
           Ta is a dictionary containing the data, quality control flag
              and attributes for the air temperature.

    Useful functions:
     pfp_utils.GetVariable(ds, label)
     pfp_utils.CreateVariable(ds, variable)

    Date: Way back in the day
    Author: PRI
    """
    def __init__(self, groups=[], global_attributes={}):
        assert isinstance(groups, list), "groups must be list"
        assert isinstance(global_attributes, dict), "global_attributes must be dict"
        self.info = {"returncodes": {"value": 0,"message": "OK"},
                      "filepath": "", "groups": ["root"],
                      "mergeserieslist": [],
                      "averageserieslist": [],
                      "intermediate": []}
        self.root = {"Attributes": copy.deepcopy(global_attributes), "Variables": {}}
        if len(groups) != 0:
            for group in groups:
                setattr(self, group, {"Attributes": {}, "Variables": {}})
                self.info["groups"].append(group)
        return