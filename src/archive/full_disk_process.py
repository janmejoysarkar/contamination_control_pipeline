#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tue Oct 21 06:21:45 PM CEST 2025
@author: sarkarjj
@hostname: hydra2

DESCRIPTION
- Runs the contamination correction on an image of choice.
- Saves the FITS file.
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import os, glob
import numpy as np
from sunpy.map import Map

def calibration(filepath, SAVE=None):
    # Open test image
    with fits.open(filepath) as test_hdu:
        test_img= test_hdu[0].data
        test_img_header= test_hdu[0].header
    FILT_NAME= test_img_header['FTR_NAME']

    # Open corresponding flat file
    if os.path.exists(flat_filename):
        with fits.open(flat_filename) as flat_hdu:
            flat= flat_hdu[0].data
    else:
        print('Check flat filepath')

    # Apply correction
    test_corr= test_img/flat
    if SAVE:
        save_fits_path=os.path.join(savepath, test_img_header['F_NAME'])
        fits.writeto(save_fits_path, test_corr, header=test_img_header, overwrite=True)
    return(test_corr)

if __name__=='__main__':
    SAVE= True # Save corrected fits?
    PLOT= True # Save preview?
    # Paths
    project_path= os.path.abspath('..')
    filt_names=['NB01', 'NB02', 'NB03', 'NB04', 'NB05', 'NB06', 'NB07', 'NB08', 'BB01', 'BB02', 'BB03']
    filt_name='NB05'
    savepath= os.path.join(project_path, 'products/full_disk')
    files= sorted(glob.glob(os.path.join(project_path, f'data/raw/*{filt_name}*'))) # Filepath for full disk images
    flat_filename= glob.glob(os.path.join(project_path, f'data/interim/flat_{files[0]}.fits'))[0]
    print(os.path.basename(flat_filename))
    for file in files:
        print(os.path.basename(file))
        corrected_image= calibration(file, SAVE)


