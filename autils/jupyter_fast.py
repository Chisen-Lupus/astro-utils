
import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, HTML
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from astropy.visualization import ImageNormalize, LinearStretch, AsinhStretch, ZScaleInterval, HistEqStretch, LogStretch, PowerDistStretch, SinhStretch, SqrtStretch, SquaredStretch
import matplotlib.patches as patches

from autils import TemporaryMatplotlibConfig, configure_inline_matplotlib, log_call

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


@log_call
def interactive_plot(data, wcs=None):
    global fig, ax, im, cbar  # Declare global variables

    def __create_triangle(x, y, l, h):
        triangle = patches.Polygon([[x, y+h], [x-l/2, y], [x+l/2, y]], 
                                closed=True, color='k')
        return triangle

    def __reset_slider(slider, value): 
        slider.value = value
    
    # with TemporaryMatplotlibConfig(backend='widget'):

    # # configure_inline_matplotlib(backend='widget')
    # print(plt.get_backend())
    # plt.imshow(np.random.rand(100, 100))
    # plt.show()

    @log_call
    def __initialize_plot(data, 
                scale_method=LinearStretch, 
                scale_range='min max', 
                scale_slider=0, 
                colormap='grey', 
                show_grid=True, 
                coordinate='world', 
                wcs=wcs):
        global fig, ax, im, cbar  # Use global variables here

        data = np.asarray(data, dtype=np.float64)
        fig = plt.figure(figsize=[8, 6])
        # parse coordinate
        if coordinate=='world':
            ax = fig.add_subplot(111, projection=wcs)
            ax.coords[0].set_axislabel_visibility_rule('ticks')
            ax.coords[0].set_ticks_visible(False)
            ax.coords[1].set_axislabel_visibility_rule('ticks')
            ax.coords[1].set_ticks_visible(False)
        elif coordinate=='image':
            ax = fig.add_subplot(111)
        else: 
            raise ValueError(f"Not a valid coordinate: {coordinate}")
        # parse scale_range
        if scale_range=='min max':
            vmin = np.min(data)
            vmax = np.max(data)
        elif scale_range=='zscale': 
            zscale_interval = ZScaleInterval()
            vmin, vmax = zscale_interval.get_limits(data)
        else: 
            raise ValueError(f'Not a valid scale_range: {scale_range}')
        # parse scale_slider
        offset = scale_slider*(vmax - vmin)
        vmax_shifted = vmax - offset
        vmin_shifted = vmin - offset
        # parse scale_method
        args = [data] if scale_method==HistEqStretch else []
        norm = ImageNormalize(stretch=scale_method(*args), 
                            vmin=vmin_shifted, vmax=vmax_shifted)
        # imshow
        im = ax.imshow(data, interpolation='none', cmap=colormap, norm=norm)
        # colorbar
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        # mark the upper and lower bound of the colorbar
        cbar.ax.add_patch(__create_triangle(0, np.min(data), 1, (vmax-vmin)/200))
        cbar.ax.add_patch(__create_triangle(1, np.min(data), 1, (vmax-vmin)/200))
        cbar.ax.add_patch(__create_triangle(0, np.max(data), 1, (vmin-vmax)/200))
        cbar.ax.add_patch(__create_triangle(1, np.max(data), 1, (vmin-vmax)/200))
        # parse show_grid
        ax.grid(show_grid)
        # show image
        fig.tight_layout()
        plt.show()
    
    @log_call
    def __setup_ax(coordinate):
        global ax
        if coordinate.value == 'world' and wcs is not None:
            ax = fig.add_subplot(111, projection=wcs)
            ax.coords[0].set_axislabel_visibility_rule('ticks')
            ax.coords[0].set_ticks_visible(False)
            ax.coords[1].set_axislabel_visibility_rule('ticks')
            ax.coords[1].set_ticks_visible(False)
        elif coordinate.value == 'image':
            ax = fig.add_subplot(111)
        else:
            raise ValueError(f"Not a valid coordinate: {coordinate.value}")

    @log_call
    def __update_coordinate(change):
        global ax, im, cbar  # Access the global variables

        # Clear current Axes to avoid creating a new plot
        ax.remove()
        # del ax
        # Update projection based on coordinate selection
        __setup_ax(coordinate)

        im = ax.imshow(data, interpolation='none')
        # Reapply plot parameters and update image
        __update_plot(None)

    @log_call
    def __update_plot(change):
        global im, cbar  # Access global variables

        # Update scale range based on selected method
        if scale_range.value == 'min max':
            vmin = np.min(data)
            vmax = np.max(data)
        elif scale_range.value == 'zscale':
            zscale_interval = ZScaleInterval()
            vmin, vmax = zscale_interval.get_limits(data)
        else:
            raise ValueError(f'Not a valid scale_range: {scale_range.value}')

        # Adjust for scale slider offset
        offset = scale_slider.value * (vmax - vmin)
        vmax_shifted = vmax - offset
        vmin_shifted = vmin - offset

        # Update normalization based on scale method
        args = [data] if scale_method.value == HistEqStretch else []
        norm = ImageNormalize(stretch=scale_method.value(*args), 
                            vmin=vmin_shifted, vmax=vmax_shifted)

        # Update image content (assuming `data` has been changed)
        im.set_data(data)  # Use `data` or `new_data` if the image data has changed
        im.set_norm(norm)  # Update the normalization
        im.set_cmap(colormap.value)  # Update colormap
        im.set_clim(vmin_shifted, vmax_shifted)  # Update color limits

        # Update colorbar to reflect changes
        cbar.update_normal(im)

        # Efficiently redraw the figure canvas
        fig.canvas.draw_idle()



    # define sliders
    
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
        value=0,
        min=-1,
        max=1,
        step=0.005,
        description='',
        continuous_update=True,
        readout=False,
        layout=Layout(width='50%')
    )

    reset_button = widgets.Button(
        description='reset',
        layout=Layout(width='10%')
    )
    reset_button.on_click(lambda b: __reset_slider(scale_slider, 0))

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

    # define ui

    ui1 = widgets.HBox([scale_method, scale_range])
    ui2 = widgets.HBox([reset_button, scale_slider, 
                        colormap, show_grid, coordinate])
    ui = widgets.VBox([ui1, ui2])

    display(ui)

    # Initialize plot

    __initialize_plot(data)

    # Link widgets to update function

    scale_method.observe(__update_plot, names='value')
    scale_range.observe(__update_plot, names='value')
    scale_slider.observe(__update_plot, names='value')
    colormap.observe(__update_plot, names='value')
    show_grid.observe(__update_plot, names='value')
    coordinate.observe(__update_coordinate, names='value')