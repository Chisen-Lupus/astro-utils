import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, HTML
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from astropy.visualization import ImageNormalize, LinearStretch, AsinhStretch, ZScaleInterval, HistEqStretch, LogStretch, PowerDistStretch, SinhStretch, SqrtStretch, SquaredStretch

np.seterr(invalid='warn')

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

display(HTML(custom_css))

class __CustomNorm(Normalize):

    def __init__(self, stretch_func=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stretch_func = stretch_func

    def __call__(self, value, clip=None):
        if self.stretch_func is not None:
            value = self.stretch_func(value)
        return super().__call__(value, clip)

def __make_plot(data, 
            scale_method, scale_range, 
            scale_slider, colormap, show_grid, coordinate, 
            wcs=None):
    data = np.asarray(data, dtype=np.float64)
    fig = plt.figure(figsize=[8, 6])
    # parse coordinate
    if coordinate=='world':
        ax = fig.add_subplot(111, projection=wcs)
        ax.coords[0].set_axislabel_visibility_rule('ticks')
        ax.coords[0].set_ticks_visible(False)
        ax.coords[1].set_axislabel_visibility_rule('ticks')
        ax.coords[1].set_ticks_visible(False)
        
    if coordinate=='image':
        ax = fig.add_subplot(111)
    else: 
        pass
    # parse scale_range
    if scale_range=='min max':
        vmin = np.min(data)
        vmax = np.max(data)
    elif scale_range=='zscale': 
        zscale_interval = ZScaleInterval()
        vmin, vmax = zscale_interval.get_limits(data)
    else: 
        pass
    # TODO: parse scale_slider
    # parse scale_method
    args = [data] if scale_method==HistEqStretch else []
    norm = ImageNormalize(stretch=scale_method(*args), vmin=vmin, vmax=vmax)


    # if coordinate
    # TODO: imshow. parse colormap here
    im = ax.imshow(data, cmap=colormap, norm=norm)
    fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    # parse show_grid
    ax.grid(show_grid)
    # show image
    fig.tight_layout()
    # plt.show()

def interactive_plot(data, wcs=None):

    scale_method = widgets.ToggleButtons(
        options=[('asinh', AsinhStretch), 
                 ('histogram', HistEqStretch), 
                 ('linear', LinearStretch), 
                 ('log', LogStretch), 
                 ('power dist', PowerDistStretch), 
                 ('sinh', SinhStretch), 
                 ('sqrt', SqrtStretch), 
                 ('squared', SquaredStretch)],
        value=LinearStretch, 
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

    show_grid = widgets.ToggleButton(
        value=True,
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
        value='world',
        layout=Layout(width='20%', height='33px'), 
        description='',
        disabled=no_wcs, 
        style={'button_width': '47%'}
    )

    ui1 = widgets.HBox([scale_method, scale_range])
    ui2 = widgets.HBox([scale_slider, colormap, show_grid, coordinate])
    ui = widgets.VBox([ui1, ui2])

    out = widgets.interactive_output(__make_plot, 
                                    {'data':  widgets.fixed(data), 
                                    'scale_method': scale_method, 
                                    'scale_range': scale_range, 
                                    'scale_slider': scale_slider, 
                                    'colormap': colormap, 
                                    'show_grid': show_grid,
                                    'coordinate': coordinate, 
                                    'wcs':  widgets.fixed(wcs)})

    display(ui, out)
    # display(out)