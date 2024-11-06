import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, HTML
import numpy as np
import matplotlib.pyplot as plt
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

class InteractivePlot:
    def __init__(self, data, wcs=None):
        self.data = np.asarray(data, dtype=np.float64)
        self.wcs = wcs
        self.vmin_shifted = None
        self.vmax_shifted = None
        self.create_widgets()
        self.initialize_plot()

    def create_widgets(self):
        # Create widgets
        self.scale_method = widgets.ToggleButtons(
            options=[('asinh', AsinhStretch), 
                     ('histogram', HistEqStretch), 
                     ('linear', LinearStretch), 
                     ('log', LogStretch), 
                     ('power dist', PowerDistStretch), 
                     ('sinh', SinhStretch), 
                     ('sqrt', SqrtStretch), 
                     ('squared', SquaredStretch)],
            value=LinearStretch, 
            layout=Layout(width='80%', height='33px')
        )

        self.scale_range = widgets.ToggleButtons(
            options=['min max', 'zscale'],
            value='min max',
            layout=Layout(width='20%', height='33px')
        )

        self.scale_slider = widgets.FloatSlider(
            value=0,
            min=-1,
            max=1,
            step=0.005,
            description='',
            continuous_update=True,
            readout=False,
            layout=Layout(width='50%')
        )

        self.colormap = widgets.Dropdown(
            options=['gist_heat', 'grey', 'seismic'],
            value='grey',
            layout=Layout(width='10%')
        )

        self.show_grid = widgets.ToggleButton(
            value=True,
            description='grid',
            disabled=False,
            button_style='',
            layout=Layout(width='10%')
        )

        self.coordinate = widgets.ToggleButtons(
            options=['image', 'world'],
            value='world',
            layout=Layout(width='20%', height='33px')
        )

        # Define layout
        ui1 = widgets.HBox([self.scale_method, self.scale_range])
        ui2 = widgets.HBox([self.scale_slider, self.colormap, self.show_grid, self.coordinate])
        display(widgets.VBox([ui1, ui2]))

        # Attach observers
        self.scale_method.observe(self.update_scale_method, names='value')
        self.scale_range.observe(self.update_scale_range, names='value')
        self.scale_slider.observe(self.update_scale_slider, names='value')
        self.colormap.observe(self.update_colormap, names='value')
        self.show_grid.observe(self.update_show_grid, names='value')
        self.coordinate.observe(self.update_coordinate, names='value')

    @log_call
    def initialize_plot(self):
        # Initial setup
        self.fig, self.ax = self.create_figure()
        self.update_axes_projection()
        self.update_image()
        self.update_colorbar()
        plt.show()

    @log_call
    def create_figure(self):
        fig = plt.figure(figsize=[8, 6])
        return fig, None

    @log_call
    def update_axes_projection(self):
        # Set up the projection on the `Axes`
        if self.coordinate.value == 'world' and self.wcs is not None:
            self.ax = self.fig.add_subplot(111, projection=self.wcs)
            self.ax.coords[0].set_axislabel_visibility_rule('ticks')
            self.ax.coords[1].set_axislabel_visibility_rule('ticks')
        else:
            self.ax = self.fig.add_subplot(111)

    @log_call
    def update_image(self):
        # Parse range and normalization
        vmin, vmax = self.calculate_vmin_vmax()
        self.vmin_shifted, self.vmax_shifted = self.apply_slider_offset(vmin, vmax)
        norm = ImageNormalize(stretch=self.scale_method.value(),
                              vmin=self.vmin_shifted, vmax=self.vmax_shifted)
        
        # Plot the image
        self.im = self.ax.imshow(self.data, cmap=self.colormap.value, norm=norm, interpolation='none')

    @log_call
    def update_colorbar(self):
        if hasattr(self, 'cbar'):
            self.cbar.remove()
        self.cbar = self.fig.colorbar(self.im, ax=self.ax, fraction=0.046, pad=0.04)

    def calculate_vmin_vmax(self):
        if self.scale_range.value == 'min max':
            return np.min(self.data), np.max(self.data)
        elif self.scale_range.value == 'zscale':
            zscale_interval = ZScaleInterval()
            return zscale_interval.get_limits(self.data)

    def apply_slider_offset(self, vmin, vmax):
        offset = self.scale_slider.value * (vmax - vmin)
        return vmin - offset, vmax - offset

    # Update methods
    def update_coordinate(self, change):
        self.ax.remove()
        self.update_axes_projection()
        self.update_image()
        self.update_colorbar()
        self.fig.canvas.draw_idle()

    def update_scale_method(self, change):
        self.update_image()
        self.update_colorbar()
        self.fig.canvas.draw_idle()

    def update_scale_range(self, change):
        self.update_image()
        self.update_colorbar()
        self.fig.canvas.draw_idle()

    def update_scale_slider(self, change):
        self.update_image()
        self.update_colorbar()
        self.fig.canvas.draw_idle()

    def update_colormap(self, change):
        self.im.set_cmap(self.colormap.value)
        self.fig.canvas.draw_idle()

    def update_show_grid(self, change):
        self.ax.grid(self.show_grid.value)
        self.fig.canvas.draw_idle()
