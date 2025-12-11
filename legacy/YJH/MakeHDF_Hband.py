# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 15:41:40 2024

@author: lvt75561
"""

from pyechelle.CCD import CCD
from pyechelle.hdfbuilder import HDFBuilder
from pyechelle.spectrograph import InteractiveZEMAX
import numpy as np
from pathlib import Path
# Open Link to a standalone OpticStudio instance
zmx = InteractiveZEMAX(name='ANDES123', zemax_filepath="HIRES_H_21jan2022_sconf.zmx")
# set basic grating specifications
zmx.set_grating(surface='ECHELLE', blaze=76)#, theta=0., gamma=0.9)
# add CCD information (only one CCD supported so far. So for instruments with multiple CCDs, you have to generate
# separate models for now.
zmx.add_ccd(1, CCD(4096, 4096, pixelsize=15))

# Add here as many fiber/fields as you wish. You don't have to fiddle with the fields in OpticStudio. The
# existing fields will be ignored/deleted.
nfiber = 75
sfi = 474 # size fiber 630 # size fiber
yfield = (np.arange(nfiber)-nfiber/2+0.5)*sfi/1000 # Vertical distribution of the fibers
xfield = np.zeros(nfiber) # Slit is vertical in its definition


for i in range(nfiber):
    # One specific fiber
    zmx.add_field(xfield[i], yfield[i], sfi, sfi, shape='circular', name='Science fiber')
    # List of diffraction order for a specific fiber
    zmx.set_orders(1, i+1, list(range(68, 83)))    

# Adjust settings for the Huygens PSF. Best to check out 'reasonable' parameters manually in ZEMAX first.
zmx.psf_settings(image_delta=3, image_sampling="128x128", pupil_sampling="64x64")

# at this point, you can interact with the spectrograph interactively, e.g. by doing something like:
#zmx.get_psf(wavelength=0.72, order=85, fiber=1, ccd_index=1)
#zmx.get_wavelength_range()
#zmx.get_transformation(wavelength=0.72, order=85, fiber=1, ccd_index=1)
# and it will/should pull the appropriate values from ZEMAX. This might be helpful for debugging.

# To generate an .HDF model file, you do:
script_dir = Path(__file__).parent
project_root = script_dir.parent
hdf_path = project_root / 'HDF' / 'ANDES_75fibre_H.hdf'
hdf = HDFBuilder(zmx, str(hdf_path))
# this will take a long time...
hdf.save_to_hdf(n_transformation_per_order=15, n_psfs_per_order=15)


hdf.close()
