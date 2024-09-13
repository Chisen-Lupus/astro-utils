import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.visualization import make_lupton_rgb, ZScaleInterval, astropy_mpl_style, AsinhStretch
from astropy.coordinates import SkyCoord
import astropy.units as u
import os

plt.style.use(astropy_mpl_style)
plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.interpolation'] = 'none'

def make_cutout(fits_file, position, radius):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        wcs = WCS(hdul[0].header)
        cutout = Cutout2D(data, position, radius, wcs=wcs, mode='partial', fill_value=None)
    return cutout

def apply_zscale(image):
    """Apply zscale normalization to an image."""
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(image)
    return np.clip((image - vmin) / (vmax - vmin), 0, 1)

def make_rgb_image(g_file, r_file, i_file, position, radius):
    g_cutout = make_cutout(g_file, position, radius)
    r_cutout = make_cutout(r_file, position, radius)
    i_cutout = make_cutout(i_file, position, radius)
    # print(r_cutout.data)
    # print(g_cutout.data.shape, r_cutout.data.shape, i_cutout.data.shape)
    
    # Apply zscale normalization to each band
    g_zscale = apply_zscale(g_cutout.data)
    r_zscale = apply_zscale(r_cutout.data)
    i_zscale = apply_zscale(i_cutout.data)
    
    # Create RGB image using the Lupton method and zscale normalized images
    rgb_image = make_lupton_rgb(i_zscale, r_zscale, g_zscale, Q=1, stretch=1)
    
    return rgb_image

def plot_grid(fits_files, coords, radius):
    # Calculate the number of rows and columns for a near-square grid
    n_images = len(fits_files)
    n_cols = int(np.ceil(np.sqrt(n_images)))
    n_rows = int(np.ceil(n_images / n_cols))
    
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 15), 
                            gridspec_kw={'wspace': 0, 'hspace': 0})  # No space between subplots

    for i, ax in enumerate(axs.flat):
        if i < len(fits_files):
            g_file, r_file, i_file = fits_files[i]
            rgb_image = make_rgb_image(g_file, r_file, i_file, coords[i], radius)
            ax.imshow(rgb_image)
            
            # Set the image border (spine color) - adds a black border to each frame
            for spine in ax.spines.values():
                spine.set_edgecolor('black')
                spine.set_linewidth(2)
            
            # Display ID as text in the lower-left corner of the image
            ax.text(0.02, 0.02, f"ID {i+1}", color='white', fontsize=15, ha='left', va='bottom', transform=ax.transAxes)
        
        # Remove axis, ticks, and grid inside the image
        ax.set_xticks([])  # Remove x-ticks
        ax.set_yticks([])  # Remove y-ticks
        ax.grid(False)     # Disable the grid
        # ax.axis('off')     # Turn off the axis lines

    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Ensure the grid fills the entire figure
    plt.show()

# # Example usage:
# fits_groups = [
#     ["id001_g.fits", "id001_r.fits", "id001_i.fits"],
#     ["id002_g.fits", "id002_r.fits", "id002_i.fits"],
#     ["id003_g.fits", "id003_r.fits", "id003_i.fits"],
#     ["id004_g.fits", "id004_r.fits", "id004_i.fits"],
#     ["id005_g.fits", "id005_r.fits", "id005_i.fits"]
# ]

# # Example coordinates for the objects
# coordinates = [
#     SkyCoord(ra=10.684*u.deg, dec=41.269*u.deg, frame='icrs'),
#     SkyCoord(ra=10.689*u.deg, dec=41.273*u.deg, frame='icrs'),
#     SkyCoord(ra=10.695*u.deg, dec=41.277*u.deg, frame='icrs'),
#     SkyCoord(ra=10.700*u.deg, dec=41.281*u.deg, frame='icrs'),
#     SkyCoord(ra=10.705*u.deg, dec=41.285*u.deg, frame='icrs')
# ]

# # Cutout radius
# cutout_radius = 10 * u.arcsec

# # Plot the grid of RGB images using zscale, borders, and IDs as text
# plot_grid(fits_groups, coordinates, cutout_radius)