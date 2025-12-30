#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Thu Oct 30 06:06:45 PM CET 2025
@author: sarkarjj
@hostname: hydra2

DESCRIPTION
"""
from sunkit_image.coalignment import calculate_match_template_shift, apply_shifts
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import astropy.units as u
import numpy as np
import os
import glob
import sunpy
import matplotlib.pyplot as plt
from astropy.io import fits

def PLOT():
    aligned_maps.peek()
    plt.show()

def get_submap(ref_map):
    """
    Get a submap from the reference image for template matching.
    Parameters:
    - ref_map: Reference SunPy map image.
    Returns:
    - Submap of the reference image.
    """
    bl= ref_map.pixel_to_world(400*u.pixel, 230*u.pixel)
    tr= ref_map.pixel_to_world(550*u.pixel, 330*u.pixel)
    ref_submap= ref_map.submap(bl, top_right=tr)
    return ref_submap
    
project_path= os.path.abspath('..')
filt_name='BB01'
files= sorted(glob.glob(os.path.join(project_path, f'data/raw/*{filt_name}*'))) # Filepath for full disk images
ref_map= sunpy.map.Map(files[0]) #First frame is taken as reference
seq = sunpy.map.Map(files[4:], sequence=True)
ref_submap= get_submap(ref_map)
o_x, o_y= [], []
ref_head=ref_submap.fits_header
ref_cdel=ref_head['CDELT1']
align_shift = calculate_match_template_shift(seq, template=ref_submap)
shift_xPix = align_shift['x'].value / ref_cdel * -1
shift_yPix = align_shift['y'].value / ref_cdel * -1
aligned_maps = apply_shifts(seq, yshift=shift_yPix * u.pixel, xshift=shift_xPix * u.pixel, clip=False)
aligned_map_arr= np.stack([m.data for m in aligned_maps], axis=0)
med= np.median(aligned_map_arr, axis=0)
flat= ref_map.data/med
