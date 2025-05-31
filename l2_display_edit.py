# Python modules
import os
from pathlib import Path
import platform
from typing import Optional
# 3rd party modules
from nicegui import events, ui
import numpy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
# PyFluxPro modules
from scripts import pfp_ck
from scripts import pfp_io
from scripts import pfp_utils

class local_file_picker(ui.dialog):

    def __init__(self, directory: str, *,
                 upper_limit: Optional[str] = ..., multiple: bool = False,
                 show_hidden_files: bool = False) -> None:
        """Local File Picker

        This is a simple file picker that allows you to select a file from the local filesystem where NiceGUI is running.

        :param directory: The directory to start in.
        :param upper_limit: The directory to stop at (None: no limit, default: same as the starting directory).
        :param multiple: Whether to allow multiple files to be selected.
        :param show_hidden_files: Whether to show hidden files.
        """
        super().__init__()

        self.path = Path(directory).expanduser()
        if upper_limit is None:
            self.upper_limit = None
        else:
            self.upper_limit = Path(directory if upper_limit == ... else upper_limit).expanduser()
        self.show_hidden_files = show_hidden_files

        with self, ui.card():
            self.add_drives_toggle()
            self.grid = ui.aggrid({
                'columnDefs': [{'field': 'name', 'headerName': 'File'}],
                'rowSelection': 'multiple' if multiple else 'single',
            }, html_columns=[0]).classes('w-96').on('cellDoubleClicked', self.handle_double_click)
            with ui.row().classes('w-full justify-end'):
                ui.button('Cancel', on_click=self.close).props('outline')
                ui.button('Ok', on_click=self._handle_ok)
        self.update_grid()

    def add_drives_toggle(self):
        if platform.system() == 'Windows':
            import win32api
            drives = win32api.GetLogicalDriveStrings().split('\000')[:-1]
            self.drives_toggle = ui.toggle(drives, value=drives[0], on_change=self.update_drive)

    def update_drive(self):
        self.path = Path(self.drives_toggle.value).expanduser()
        self.update_grid()

    def update_grid(self) -> None:
        paths = list(self.path.glob('*'))
        if not self.show_hidden_files:
            paths = [p for p in paths if not p.name.startswith('.')]
        paths.sort(key=lambda p: p.name.lower())
        paths.sort(key=lambda p: not p.is_dir())

        self.grid.options['rowData'] = [
            {
                'name': f'📁 <strong>{p.name}</strong>' if p.is_dir() else p.name,
                'path': str(p),
            }
            for p in paths
        ]
        if (self.upper_limit is None and self.path != self.path.parent) or \
                (self.upper_limit is not None and self.path != self.upper_limit):
            self.grid.options['rowData'].insert(0, {
                'name': '📁 <strong>..</strong>',
                'path': str(self.path.parent),
            })
        self.grid.update()

    def handle_double_click(self, e: events.GenericEventArguments) -> None:
        self.path = Path(e.args['data']['path'])
        if self.path.is_dir():
            self.update_grid()
        else:
            self.submit([str(self.path)])

    async def _handle_ok(self):
        rows = await self.grid.get_selected_rows()
        self.submit([r['path'] for r in rows])

async def main_file_select():
    """
    Purpose:
     Waits until the user selects a control file using local_file_picker()
     and then reads the input and output netCDF files before handing off
     to the display_data() function.
    Called by: main_gui()
    Calls: display_data()
    Date: May 2025
    Author: PRI
    """
    # code needs to be refactored into classes to avoid use of global
    global cfg, ds_l1, ds_l2
    result = await local_file_picker('~', multiple=True)
    if result:
        ui.notify(f"Selected file: " + result[0])
        cfg = pfp_io.get_controlfilecontents(result[0], mode="verbose")
        data_path = cfg["Files"]["file_path"]
        in_filename = cfg["Files"]["in_filename"]
        out_filename = cfg["Files"]["out_filename"]
        in_uri = os.path.join(data_path, in_filename)
        out_uri = os.path.join(data_path, out_filename)
        ds_l1 = pfp_io.NetCDFRead(in_uri)
        ds_l2 = pfp_io.NetCDFRead(out_uri)
        display_data()

def display_data():
    """
    Purpose:
     Displays a list of variables and variable groups as tree structure
     on the left of the browser window and waits for the user to select
     an individual variable or a group of variables for visualisation.
    Called by: main_file_select()
    Calls: on_select() when the user chooses a group or a variable.
    Date: May 2025
    Author: PRI
    """
    labels = sorted(list(cfg["Variables"].keys()))
    for label in list(labels):
        if "DependencyCheck" in cfg["Variables"][label]:
            source = cfg["Variables"][label]["DependencyCheck"]["source"]
            source = source.split(",")
            for src in source:
                if src in labels:
                    labels.remove(src)
    groups = {"Groups": ["Radiation", "Fluxes", "Meteorology", "Covariances",
                         "Soil temperature", "Soil moisture", "Soil heat flux"],
              "Variables": labels}
    tree_data = [dict_to_tree_data(groups, label="Data")]
    tree_card.clear()
    with tree_card:
        ui.tree(tree_data, label_key='id', on_select=lambda e: on_select(e))
    return

def on_select(e):
    """
    Purpose:
     Checks to see what the user has selected in the tree structure and then
     perorms the approrpiate action.
    Called by: display_data()
    Calls:
    Date: May 2025
    Author: PRI
    """
    ui.notify(e.value)
    ds_l2_labels = list(ds_l2.root["Variables"].keys())
    plt_labels = []
    # did the user select an item from the Group section of the tree?
    if e.value in ["Radiation", "Fluxes", "Meteorology", "Covariances",
                   "Soil temperature", "Soil moisture", "Soil heat flux"]:
        # if so, what group was selected?
        if e.value == "Radiation":
            # get a list of radiation variables
            for item in ["Fsd", "Fsu", "Fld", "Flu", "Fn"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                # remove any radiation variables that appear as a dependency
                # not sure why I did this now ...
                for label in list(l):
                    if "DependencyCheck" in cfg["Variables"][label]:
                        source = cfg["Variables"][label]["DependencyCheck"]["source"]
                        source = source.split(",")
                        for src in source:
                            if src in l:
                                l.remove(src)
                plt_labels = plt_labels + l
        elif e.value == "Fluxes":
            # get a list of flux variables
            for item in ["Fco2", "Fe", "Fh", "Fm"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                for label in list(l):
                    if "DependencyCheck" in cfg["Variables"][label]:
                        source = cfg["Variables"][label]["DependencyCheck"]["source"]
                        source = source.split(",")
                        for src in source:
                            if src in l:
                                l.remove(src)
                plt_labels = plt_labels + l
        elif e.value == "Meteorology":
            # get a list of meteorological variables
            for item in ["Ta", "AH", "RH", "Ws", "Wd"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                plt_labels = plt_labels + l
        elif e.value == "Covariances":
            # get a list of covariance variables
            for item in ["UxA", "UxC", "UxT", "UyA", "UyC", "UyT", "UzA", "UzC", "UzT"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                plt_labels = plt_labels + l
        elif e.value == "Soil temperature":
            # get a list of soil temperature variables
            for item in ["Ts"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                plt_labels = plt_labels + l
        elif e.value == "Soil moisture":
            # get a list of soil moisture variables
            for item in ["Sws"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                plt_labels = plt_labels + l
        elif e.value == "Soil heat flux":
            # get a list of soil heat filled variables
            for item in ["Fg"]:
                l = [i for i in ds_l2_labels if i[:len(item)] == item and "_QCFlag" not in i]
                plt_labels = plt_labels + l
        # now plot a stacked time series of variables in the selected group
        with plot_card:
            if len(plt_labels) > 0:
                plot_card.clear()
                # get a row for the plot
                container = ui.row()
                with container:
                    # put up buttons to allow the user to select between stacked time series or
                    # fingerprint plots
                    ui.button(text='Time series',
                              on_click=lambda e, p=plt_labels, c=container: onclick_plot_time_series_stacked(e, p, c),
                              color='primary')
                    ui.button(text='Fingerprint',
                              on_click=lambda e, p=plt_labels, c=container: onclick_plot_fingerprint(e, p, c),
                              color='primary')
                    # and plot teh stacked time series by default
                    plot_time_series_stacked(plt_labels)
            else:
                ui.notify('Nothing to plot (' + e.value + ')')
    # did the user select a single variable from the tree?
    elif e.value in ds_l2_labels:
        # if so, plot it
        plt_label = e.value
        vcfg = cfg["Variables"][plt_label]
        with plot_card:
            plot_card.clear()
            # use different plot if this variable has dependencies
            if "DependencyCheck" in vcfg:
                dependencies = vcfg["DependencyCheck"]["source"].split(',')
                plt_labels = [plt_label] + dependencies
                plot_time_series_dependency(plt_labels)
            else:
                plot_time_series_single(plt_label)
    return
def onclick_plot_time_series_stacked(e, p, c):
    c.remove(2)
    plot_time_series_stacked(p)
    return
def onclick_plot_fingerprint(e, p, c):
    c.remove(2)
    plot_fingerprint(p)
    return
def dict_to_tree_data(data, label="Data"):
    """Converts a Python dictionary (or list) to the list format needed for ui.tree."""
    if isinstance(data, dict):
        children = [dict_to_tree_data(value, key) for key, value in data.items()]
        return {"id": label, "children": children}
    elif isinstance(data, list):
        children = [{"id": str(item)} for item in data]
        return {"id": label, "children": children}
    else:  # Handle other data types (e.g., strings, numbers)
        return {"id": f"{label}: {data}"}
def plot_fingerprint(labels):
    """
    Purpose:
     Plot fingerprints of the selected variables.
     Mostly stolen from PyFluxPro.
    Called by: onclick_plot_fingerprint()
    Calls: pfp_utils.GetVariable()
           plotly.subplots.make_subplots()
    Date: May 2025
    Author: PRI
    """
    ncols = len(labels)
    nrows = 1
    ts = int(ds_l2.root["Attributes"]["time_step"])
    nperday = int(24/(ts/60))
    dt = pfp_utils.GetVariable(ds_l2, "DateTime", match="wholedays")
    ndays = int(len(dt["Data"])/nperday)
    data = {}
    for l in labels:
        var = pfp_utils.GetVariable(ds_l2, l, match="wholedays")
        data[l] = numpy.ma.filled(var["Data"], fill_value=numpy.nan)
    fig = make_subplots(rows=nrows, cols=ncols, shared_yaxes=True)
    for n, l in enumerate(labels):
        fig.add_trace(go.Heatmap(z=data[l].reshape(ndays, nperday),
                                 colorscale="jet", showscale=False),
                      row=1, col=n+1)
    fig.update_layout(height=700, width=900, autosize=False)
    return ui.plotly(fig)
def plot_time_series_stacked(labels):
    """
    Purpose:
     Plot stacked time series of the selected variables.
     Mostly stolen from PyFluxPro.
    Called by: onclick_plot_time_series_stacked()
    Calls: pfp_utils.GetVariable()
           plotly.subplots.make_subplots()
    Date: May 2025
    Author: PRI
    """
    nrows = len(labels)
    ncols = 1
    fig = make_subplots(rows=nrows, cols=ncols, shared_xaxes=True, vertical_spacing=0.02)
    for n, label in enumerate(labels):
        var = pfp_utils.GetVariable(ds_l2, label)
        var["DateTime"] = var["DateTime"].astype('datetime64[ns]')
        var["Data"] = numpy.ma.filled(var["Data"], fill_value=numpy.nan)
        fig.add_trace(go.Scatter(x=var["DateTime"], y=var["Data"], connectgaps=False),
                      row=n+1, col=1)
    site_name = ds_l2.root["Attributes"]["site_name"]
    start_date = ds_l2.root["Attributes"]["time_coverage_start"]
    end_date = ds_l2.root["Attributes"]["time_coverage_end"]
    title = site_name + ": " + start_date + " to " + end_date
    fig.update_layout(height=700, width=900, title_text=title, showlegend=False)
    return ui.plotly(fig)
@ui.refreshable
def plot_time_series_dependency(labels):
    """
    Purpose:
     Plot time series of the selected variable with dependencies.
     Mostly stolen from PyFluxPro.
    Called by: on_select()
    Calls: pfp_utils.GetVariable()
           plotly.subplots.make_subplots()
    Date: May 2025
    Author: PRI
    """
    l0 = labels[0]
    vcfg = cfg["Variables"][l0]
    nrows = len(labels)
    ncols = 1
    fig = make_subplots(rows=nrows, cols=ncols, shared_xaxes=True, vertical_spacing=0.02)
    for n, label in enumerate(labels):
        var_in = pfp_utils.GetVariable(ds_l1, label)
        var_in["DateTime"] = var_in["DateTime"].astype('datetime64[ns]')
        var_in["Data"] = numpy.ma.filled(var_in["Data"], fill_value=numpy.nan)
        var_out = pfp_utils.GetVariable(ds_l2, label)
        var_out["DateTime"] = var_out["DateTime"].astype('datetime64[ns]')
        var_out["Data"] = numpy.ma.filled(var_out["Data"], fill_value=numpy.nan)
        line_in = {"color": ('rgba(255,0,0,0.5)')}
        fig.add_trace(go.Scatter(x=var_in["DateTime"], y=var_in["Data"], connectgaps=False,
                                 line=line_in), row=n+1, col=1)
        line_out = {"color": ('rgba(0,0,255,1)')}
        fig.add_trace(go.Scatter(x=var_out["DateTime"], y=var_out["Data"], connectgaps=False,
                                 line=line_out), row=n+1, col=1)
    site_name = ds_l2.root["Attributes"]["site_name"]
    start_date = ds_l2.root["Attributes"]["time_coverage_start"]
    end_date = ds_l2.root["Attributes"]["time_coverage_end"]
    title = site_name + ": " + start_date + " to " + end_date
    fig.update_layout(height=700, width=800, title_text=title, showlegend=False)
    with ui.row():
        ui.plotly(fig)
        with ui.column():
            ui.button('Apply', on_click=lambda e, l=l0: refresh_qc_plot(e, l))
            if "RangeCheck" in vcfg:
                lwr = float(vcfg["RangeCheck"]["lower"])
                upr = float(vcfg["RangeCheck"]["upper"])
                with ui.card():
                    ui.label(l0).style('width: 100px')
                    with ui.row():
                        ui.number(label='Lower', value=lwr, on_change=lambda e, l=l0: lwr_changed(e, l)).style('width: 50px')
                        ui.number(label='Upper', value=upr, on_change=lambda e, l=l0: upr_changed(e, l)).style('width: 50px')
            if "DependencyCheck" in vcfg:
                source = str(vcfg["DependencyCheck"]["source"])
                for d in source.split(","):
                    dcfg = cfg["Variables"][d]
                    dlwr = float(dcfg["RangeCheck"]["lower"])
                    dupr = float(dcfg["RangeCheck"]["upper"])
                    with ui.card():
                        ui.label(d).style('width: 100px')
                        with ui.row():
                            ui.number(label='Lower', value=dlwr, on_change=lambda e, l=d: lwr_changed(e, l)).style('width: 50px')
                            ui.number(label='Upper', value=dupr, on_change=lambda e, l=d: upr_changed(e, l)).style('width: 50px')
    return
@ui.refreshable
def plot_time_series_single(label):
    vcfg = cfg["Variables"][label]
    var_l1 = pfp_utils.GetVariable(ds_l1, label)
    var_l1["DateTime"] = var_l1["DateTime"].astype('datetime64[ns]')
    var_l1["Data"] = numpy.ma.filled(var_l1["Data"], fill_value=numpy.nan)
    var_l2 = pfp_utils.GetVariable(ds_l2, label)
    var_l2["DateTime"] = var_l2["DateTime"].astype('datetime64[ns]')
    var_l2["Data"] = numpy.ma.filled(var_l2["Data"], fill_value=numpy.nan)

    line_in = {"color": ('rgba(255,0,0,0.5)')}
    fig = make_subplots(rows=1, cols=1, shared_xaxes=True, vertical_spacing=0.02,
                        specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=var_l1["DateTime"], y=var_l1["Data"], connectgaps=False,
                             line=line_in, name="in"), secondary_y=False, row=1, col=1)
    line_out = {"color": ('rgba(0,0,255,1)')}
    fig.add_trace(go.Scatter(x=var_l2["DateTime"], y=var_l2["Data"], connectgaps=False,
                             line=line_out, name="out"), secondary_y=True, row=1, col=1)
    site_name = ds_l2.root["Attributes"]["site_name"]
    start_date = ds_l2.root["Attributes"]["time_coverage_start"]
    end_date = ds_l2.root["Attributes"]["time_coverage_end"]
    title = site_name + ": " + start_date + " to " + end_date
    fig.update_layout(height=750, width=800, title_text=title, showlegend=False)
    with ui.row():
        ui.plotly(fig)
        with ui.column():
            ui.button('Apply', on_click=lambda e, l=label: refresh_qc_plot(e, l))
            if "RangeCheck" in vcfg:
                lwr = float(vcfg["RangeCheck"]["lower"])
                upr = float(vcfg["RangeCheck"]["upper"])
                with ui.row():
                    ui.number(label='Lower', value=lwr, on_change=lambda e, l=label: lwr_changed(e, l)).style('width: 50px')
                    ui.number(label='Upper', value=upr, on_change=lambda e, l=label: upr_changed(e, l)).style('width: 50px')
            else:
                with ui.row():
                    ui.number(label='Lower', on_change=lambda e, l=label: lwr_changed(e, l)).style('width: 50px')
                    ui.number(label='Upper', on_change=lambda e, l=label: upr_changed(e, l)).style('width: 50px')
            with ui.row():
                with ui.input("From", on_change=lambda e, l=label: from_changed(e, l)).style('width: 125px') as from_date:
                    with ui.menu().props('no-parent-event') as menu:
                        with ui.date().bind_value(from_date):
                            with ui.row().classes('justify-end'):
                                ui.button('Close', on_click=menu.close).props('flat')
                    with from_date.add_slot('append'):
                        ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')
            with ui.row():
                with ui.input("To", on_change=lambda e, l=label: to_changed(e, l)).style('width: 125px') as to_date:
                    with ui.menu().props('no-parent-event') as menu:
                        with ui.date().bind_value(to_date):
                            with ui.row().classes('justify-end'):
                                ui.button('Close', on_click=menu.close).props('flat')
                    with to_date.add_slot('append'):
                        ui.icon('edit_calendar').on('click', menu.open).classes('cursor-pointer')
            #if "ExcludeDates" in vcfg:
                #with ui.row():
                    #with ui.scroll_area().classes('w-32 h-32 border'):
                        #for key in list(vcfg["ExcludeDates"].keys()):
                            #l = str(key) + " " + vcfg["ExcludeDates"][key]
                            #ui.label(l)
    return
def lwr_changed(e, l):
    cfg["Variables"][l]["RangeCheck"]["lower"] = e.value
    return
def upr_changed(e, l):
    cfg["Variables"][l]["RangeCheck"]["upper"] = e.value
    return
def from_changed(e, l):
    if "ExcludeDates" not in cfg["Variables"][l]:
        cfg["Variables"][l]["ExcludeDates"] = {}
    n = len(cfg["Variables"][l]["ExcludeDates"].keys())
    cfg["Variables"][l]["ExcludeDates"][str(n)] = e.value + " 00:30"
    return
def to_changed(e, l):
    if "ExcludeDates" not in cfg["Variables"][l]:
        cfg["Variables"][l]["ExcludeDates"] = {}
    n = len(cfg["Variables"][l]["ExcludeDates"].keys())
    cfg["Variables"][l]["ExcludeDates"][str(n - 1)] = cfg["Variables"][l]["ExcludeDates"][str(n - 1)] + "," + e.value + " 00:00"
    return
def refresh_qc_plot(e, l):
    vcfg = cfg["Variables"][l]
    if "RangeCheck" in vcfg:
        pfp_ck.do_rangecheck(cfg, {"in": ds_l1, "out": ds_l2}, "Variables", l)
    if "ExcludeDates" in vcfg:
        pfp_ck.do_excludedates(cfg, {"in": ds_l1, "out": ds_l2}, "Variables", l)
    if "DependencyCheck" in vcfg:
        pfp_ck.do_dependencycheck(cfg, {"in": ds_l1, "out": ds_l2}, "Variables", l)
        plot_time_series_dependency.refresh()
    else:
        plot_time_series_single.refresh()
    return
def run_current():
    pass
def save_file():
    pass
def main_gui():
    with ui.row():
        ui.button('Choose file', on_click=main_file_select, icon='folder')
        ui.button('Run', on_click=run_current, icon='folder')
        ui.button('Save file', on_click=save_file, icon='folder')
    with ui.row():
        global tree_card, plot_card
        tree_card = ui.card().style('height: 800px; width: 250px; overflow: auto;')
        plot_card = ui.card().style('height: 800px; width: 1000px;')
    return

if __name__ in {"__main__", "__mp_main__"}:
    @ui.page('/')
    def index():
        main_gui()
    ui.run()
    #ui.run(on_air='K1Wegv0E2UurtNPF', reconnect_timeout=10)