from l2_display_edit import main_gui
#from nicegui_plotly_test import main
from nicegui import ui

@ui.page('/')
def main_page() -> None:
    main_gui()

ui.run(show=False, port=8050)