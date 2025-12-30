#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tue Nov 18 03:04:32 PM CET 2025
@author: janmejoy
@hostname: machine

DESCRIPTION

- Made to align and cut SUIT images to see the same frame as IRIS SJI convolved images.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob, os
from sunpy.map import Map
import astropy.units as u
from sunkit_image.coalignment import apply_shifts, calculate_match_template_shift

def normalize_image(img_data):
    return (img_data-np.nanmin(img_data))/(np.nanmax(img_data)-np.nanmin(img_data))

def visualize(im0, im1):
    """
    DESCRIPTION:
    See a preview of the individual images.
    INPUT:
    - map1, map3: Sunpy Map
    - flatframe: numpy 2D array
    RETURNS: Visualization.
    """
    VMN= 0
    VMX= 1
    fig, ax= plt.subplots(1,2, sharex=True, sharey=True)
    im0= ax[0].imshow(im0, origin='lower', vmin=VMN, vmax=VMX)
    plt.colorbar(im0, ax=ax[0])
    im1= ax[1].imshow(im1, origin='lower', vmin=VMN, vmax=VMX)
    plt.colorbar(im1, ax=ax[1])
    plt.show()

def center_coordinate(bl, tr):
    bl_arcsec = (bl.Tx.to('arcsec').value, bl.Ty.to('arcsec').value)
    tr_arcsec = (tr.Tx.to('arcsec').value, tr.Ty.to('arcsec').value)
    x= np.round((bl_arcsec[0]+tr_arcsec[0])/2, 1)
    y= np.round((bl_arcsec[1]+tr_arcsec[1])/2, 1)
    print(f'Center coordinates (X,Y): {x}", {y}"')

if __name__=='__main__':
    PLOT= True
    SAVE= True 
    project_path= os.path.abspath("..")
    plot_savepath= os.path.join(project_path, 'reports/validation_with_iris/')
    iris_filename= os.path.join(project_path, 'data/external/2025_06_02_SJI.fits')
    suit_filename= os.path.join(project_path, 'data/processed/SUT_T25_0800_001000_Lev1.0_2025-06-02T05.41.03.860_0973NB03.fits')

    iris_map= Map(iris_filename)
    suit_map= Map(suit_filename)
    suit_map= suit_map.rotate(suit_map.meta['CROTA2']*u.deg)

    bl= iris_map.bottom_left_coord
    tr= iris_map.top_right_coord

    bl_suit= bl.transform_to(suit_map.coordinate_frame)
    tr_suit= tr.transform_to(suit_map.coordinate_frame)
    suit_submap= suit_map.submap(bottom_left=bl_suit, top_right=tr_suit)

    iris_map= Map(iris_map.data[:, 15:-15], iris_map.meta)
    suit_submap= Map(suit_submap.data[:, 15:-15], suit_submap.meta)

    iris_normalized= normalize_image(iris_map.data)
    suit_normalized= normalize_image(suit_submap.data)
    center_coordinate(bl, tr)
    ##################

    if PLOT:
        ## Histogram
        plt.figure(figsize=(12,5))
        plt.hist(np.ravel(iris_normalized), bins=1000, label="IRIS")
        plt.hist(np.ravel(suit_normalized), bins=1000, alpha=0.8, label="SUIT")
        plt.legend()
        plt.xlabel('Relative intensity')
        plt.ylabel('No. of pixels')
        plt.title(f'{os.path.basename(iris_filename)} and {os.path.basename(suit_filename)[27:-5]}')
        if SAVE: plt.savefig(os.path.join(plot_savepath, 'histogram.png'), dpi=300)
        plt.show()

        ## Rel Intensity Plot
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.imshow(iris_normalized, origin='lower')
        plt.title(os.path.basename(iris_filename))
        plt.colorbar()
        plt.subplot(1,2,2)
        plt.imshow(suit_normalized, origin='lower')
        plt.title(os.path.basename(suit_filename)[27:-5])
        plt.colorbar()
        if SAVE: plt.savefig(os.path.join(plot_savepath, 'rel_intensity_plot.png'), dpi=300)
        plt.show()

        ## Sunpy map
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1, projection=iris_map)
        ax2 = fig.add_subplot(1, 2, 2, projection=suit_map)
        iris_map.plot(axes=ax1)
        suit_submap.plot(axes=ax2)
        plt.show()

        ## Scatter plot
        suit_data= np.ravel(suit_normalized) 
        iris_data= np.ravel(iris_normalized)[:suit_data.shape[0]]
        colors= np.arange(iris_data.shape[0])
        plt.figure()
        plt.hist2d(iris_data, suit_data, bins=100)
        plt.show()
