import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, HTML
import numpy as np
import matplotlib.pyplot as plt

def __make_plot(data, 
            scale_method, scale_range, 
            scale_slider, colormap, grid, coordinate, 
            wcs=None):
    plt.figure(figsize=[8, 6])
    # parse scaling method
    if scale_method=='linear': 
        strech = lambda x: x
    elif scale_method=='log':
        strech = np.log
    elif scale_method=='power':
        strech = np.exp
    elif scale_method=='sqrt':
        strech = np.sqrt
    elif scale_method=='squared':
        strech = lambda x: x**2
    elif scale_method=='asinh': 
        strech = np.arcsinh
    elif scale_method=='sinh': 
        strech = np.sinh
    elif scale_method=='histogram': 
        strech = NotImplemented
    else: 
        strech = NotImplemented
    # parse scale_range and scale_slider
    vmin = NotImplemented
    vmax = NotImplemented
    # TODO: parse coordinate
    # TODO: parse wcs
    has_wcs = wcs is not None
    # TODO: temp imshow
    plt.imshow(strech(data)*scale_slider*2, vmin=np.min(strech(data)), vmax=np.max(strech(data)), cmap=colormap)
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    # parse grid
    if grid: 
        plt.grid('off')
    else: 
        plt.grid()
    plt.tight_layout()
    plt.show()


def interactive_plot(data, wcs=None):

    scale_method = widgets.ToggleButtons(
        options=['linear', 'log', 'power', 'sqrt', 'squared', 'asinh', 'sinh', 'histogram'],
        value='linear', 
        layout=Layout(width='80%', height='33px'), 
        description='',
        disabled=False, 
        style={'button_width': '11.8%'}
    )

    scale_range = widgets.ToggleButtons(
        options=['min max', 'zscale'],
        value='min max',
        layout=Layout(width='20%', height='33px'), 
        description='',
        disabled=False, 
        style={'button_width': '47%'}
    )

    scale_slider = widgets.FloatSlider(
        value=0.5,
        min=0,
        max=1,
        step=0.005,
        description='',
        continuous_update=True,
        readout=False,
        layout=Layout(width='60%')
    )

    colormap = widgets.Dropdown(
        options=['gist_heat', 'grey', 'seismic'],
        value='grey',
        description='',
        disabled=False,
        layout=Layout(width='10%'), 
    )

    grid = widgets.ToggleButton(
        value=False,
        description='grid',
        disabled=False,
        button_style='',
        tooltip='Description',
        icon='check',
        layout=Layout(width='10%'), 
    )

    no_wcs = wcs is None
    coordinate = widgets.ToggleButtons(
        options=['image', 'world'],
        value='image',
        layout=Layout(width='20.2%', height='33px'), 
        description='',
        disabled=no_wcs, 
        style={'button_width': '47%'}
    )

    custom_css = """
    <style>
        .widget-slider .slider {
            height: 16px; /* Thickness of the slider bar */
        }
        .widget-slider .noUi-handle {
            margin-top: 6px; /* Adjust to center the handle vertically */
        }
    </style>
    """

    ui1 = widgets.HBox([scale_method, scale_range])
    ui2 = widgets.HBox([scale_slider, colormap, grid, coordinate])
    ui = widgets.VBox([ui1, ui2])

    out = widgets.interactive_output(__make_plot, 
                                    {'data':  widgets.fixed(data), 
                                    'scale_method': scale_method, 
                                    'scale_range': scale_range, 
                                    'scale_slider': scale_slider, 
                                    'colormap': colormap, 
                                    'grid': grid,
                                    'coordinate': coordinate, 
                                    'wcs':  widgets.fixed(wcs)})

    display(HTML(custom_css), ui, out)