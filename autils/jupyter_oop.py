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
    /* Thickness of the slider bar */
    .widget-slider .slider {
        height: 16px; 
    }
    /* Adjust to center the handle vertically */
    .widget-slider .noUi-handle {
        margin-top: 6px; 
    }
</style>
"""

display(HTML(custom_css))

class InteractivePlot:
    def __init__(self, data, wcs=None):
        # parse data
        self.data = np.asarray(data, dtype=np.float64)
        self.wcs = wcs
        # create widgets
        self.scale_method = None
        self.scale_range = None
        self.scale_slider = None
        self.reset_button = None
        self.reset_button = None
        self.colormap = None
        self.show_grid = None
        self.coordinate = None
        self.__create_widgets()
        # main plotting objects
        self.fig = None
        self.ax = None
        self.im = None
        self.cbar = None
        # helper plotting parameters
        self.norm = None
        self.vmin = None
        self.vmax = None
        self.vmin_shifted = None
        self.vmax_shifted = None
        # plot as soon as the class is initialized
        self.initialize_plot()

    def __create_widgets(self):
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
            layout=Layout(width='80%', height='33px'), 
            description='',
            disabled=False, 
            style={'button_width': '11.8%'}
        )

        self.scale_range = widgets.ToggleButtons(
            options=['min max', 'zscale'],
            value='min max',
            layout=Layout(width='20%', height='33px'), 
            description='',
            disabled=False, 
            style={'button_width': '47%'}
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

        self.reset_button = widgets.Button(
            description='reset',
            layout=Layout(width='10%')
        )
        self.reset_button.on_click(lambda b: self.__reset_slider(self.scale_slider, 0))

        self.colormap = widgets.Dropdown(
            options=['gist_heat', 'grey', 'seismic'],
            value='grey',
            description='',
            disabled=False,
            layout=Layout(width='10%'), 
        )

        self.show_grid = widgets.ToggleButton(
            value=True,
            description='grid',
            disabled=False,
            button_style='',
            tooltip='Description',
            icon='check',
            layout=Layout(width='10%'), 
        )

        self.coordinate = widgets.ToggleButtons(
            options=['image', 'world'],
            value='world',
            layout=Layout(width='20%', height='33px'), 
            description='',
            disabled=self.wcs is None, 
            style={'button_width': '47%'}
        )

        ui1 = widgets.HBox([self.scale_method, self.scale_range])
        ui2 = widgets.HBox([self.reset_button, self.scale_slider, 
                            self.colormap, self.show_grid, self.coordinate])
        ui = widgets.VBox([ui1, ui2])

        display(ui)

        # Attach observers
        self.scale_method.observe(self.update_scale_method, names='value')
        self.scale_range.observe(self.update_scale_range, names='value')
        self.scale_slider.observe(self.update_scale_slider, names='value')
        self.colormap.observe(self.update_colormap, names='value')
        self.show_grid.observe(self.update_show_grid, names='value')
        self.coordinate.observe(self.update_coordinate, names='value')

    @log_call
    def initialize_plot(self):
        self.__setup_fig()
        # parse coordinate
        self.__setup_ax()
        # parse scale_range
        if self.scale_range.value=='min max':
            self.vmin = np.min(self.data)
            self.vmax = np.max(self.data)
        elif self.scale_range.value=='zscale': 
            zscale_interval = ZScaleInterval()
            self.vmin, self.vmax = zscale_interval.get_limits(self.data)
        else: 
            raise ValueError(f'Not a valid scale_range: {self.scale_range}')
        # parse scale_slider
        offset = self.scale_slider.value*(self.vmax - self.vmin)
        self.vmax_shifted = self.vmax - offset
        self.vmin_shifted = self.vmin - offset
        # parse scale_method
        args = [self.data] if self.scale_method.value==HistEqStretch else []
        self.norm = ImageNormalize(stretch=self.scale_method.value(*args), 
                                   vmin=self.vmin_shifted, vmax=self.vmax_shifted)
        # imshow
        self.__setup_im()
        # colorbar
        self.__setup_cbar()
        # parse show_grid
        self.ax.grid(self.show_grid.value)
        # show image
        self.fig.tight_layout()
        plt.show()

    @log_call
    def __setup_fig(self):
        self.fig = plt.figure(figsize=[8, 6])

    @log_call
    def __setup_ax(self):
        if self.coordinate.value=='world' and self.wcs is not None:
            self.ax = self.fig.add_subplot(111, projection=self.wcs)
            self.ax.coords[0].set_axislabel_visibility_rule('ticks')
            self.ax.coords[0].set_ticks_visible(False)
            self.ax.coords[1].set_axislabel_visibility_rule('ticks')
            self.ax.coords[1].set_ticks_visible(False)
        elif self.coordinate.value=='image':
            self.ax = self.fig.add_subplot(111)
        else:
            raise ValueError(f"Not a valid coordinate: {self.coordinate.value}")

    @log_call
    def __setup_im(self):
        self.im = self.ax.imshow(self.data, 
                                 interpolation='none', 
                                 cmap=self.colormap.value, 
                                 norm=self.norm)

    @log_call
    def __setup_cbar(self):
        self.cbar = self.fig.colorbar(self.im, ax=self.ax, fraction=0.046, pad=0.04)
        # mark the upper and lower bound of the colorbar
        self.cbar.ax.add_patch(self.__create_triangle(0, np.min(self.data), 
                                                      1, (self.vmax-self.vmin)/200))
        self.cbar.ax.add_patch(self.__create_triangle(1, np.min(self.data), 
                                                      1, (self.vmax-self.vmin)/200))
        self.cbar.ax.add_patch(self.__create_triangle(0, np.max(self.data), 
                                                      1, (self.vmin-self.vmax)/200))
        self.cbar.ax.add_patch(self.__create_triangle(1, np.max(self.data), 
                                                      1, (self.vmin-self.vmax)/200))










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

    # helper functions

    def __create_triangle(self, x, y, l, h):
        triangle = patches.Polygon([[x, y+h], [x-l/2, y], [x+l/2, y]], 
                                closed=True, color='k')
        return triangle

    def __reset_slider(self, slider, value): 
        slider.value = value