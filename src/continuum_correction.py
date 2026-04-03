#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2025-12-30 12:13:27
@author: janmejoyarch
@hostname: suitpoc1

DESCRIPTION
SUIT pipeline implementation of contamination correction.
Made for continuum channel images only.
"""
import os
import glob
import numpy as np
import multiprocessing
import astropy.units as u
from sunpy.map import Map
from astropy.io import fits
from datetime import datetime, UTC
from sunkit_image.coalignment import apply_shifts
from astropy.convolution import convolve, Box2DKernel
from sunpy.map.maputils import all_coordinates_from_map, coordinate_is_on_solar_disk

def blur(data, kernel): #blurring function
    return(convolve(data, Box2DKernel(kernel), normalize_kernel=True))

def alignmaps(files):
    '''
    Makes flat frame from a series of images
    '''
    m= Map(files[0])
    if ('enable' in m.meta['BIN_EN']):
        files= files[:10]
    seq = Map(files, sequence=True)
    template_map= seq[0]
    x_arry=[template_map.meta['CRPIX1']-m.meta['CRPIX1'] for i, m in enumerate(seq)]
    y_arry=[template_map.meta['CRPIX2']-m.meta['CRPIX2'] for i, m in enumerate(seq)]
    aligned_maps = apply_shifts(seq, yshift= y_arry * u.pixel, xshift=x_arry * u.pixel, clip=False)
    return aligned_maps, seq

def makeflat(aligned_maps, seq):
    map_arr= np.stack([m.data for m in seq], axis=0)
    raw_med= np.median(map_arr, axis=0)
    aligned_imgs= np.stack([m.data for m in aligned_maps], axis=0)
    med= np.median(aligned_imgs, axis=0)
    med[med==0]=1
    flat_frame= raw_med/med
    flat_frame=flat_frame/blur(flat_frame, 25) # High pass filtering
    flat_frame[flat_frame==0]=1
    header=fits.Header()
    header['NAXIS1']= aligned_maps[0].meta['NAXIS1']
    header['NAXIS2']= aligned_maps[0].meta['NAXIS2']
    header['FTR_NAME']= aligned_maps[0].meta['FTR_NAME']
    header['T_OBS']= aligned_maps[0].meta['T_OBS']
    header['WAVELNTH']= aligned_maps[0].meta['WAVELNTH']
    header['GRT_DT']= str(datetime.now(UTC))
    if SAVE:
        flat_savepath= os.path.join(project_path, f'data/interim/flat.fits')
        fits.writeto(flat_savepath, flat_frame, header, overwrite=True)
    return flat_frame

def fd_correction(file):
    '''
    Applies correction using available image
    '''
    m= Map(file)
    hpc_coords= all_coordinates_from_map(m)
    mask= np.invert(coordinate_is_on_solar_disk(hpc_coords))
    flat_frame[mask]=1
    corrected_img_data= m.data/flat_frame
    corrected_img_data= np.nan_to_num(corrected_img_data, nan=0.0)
    corrected_img_data[corrected_img_data> 6e4]=0
    corrected_map= Map(corrected_img_data, m.meta)
    if SAVE:
        filename= m.meta['F_NAME']
        img_savepath= os.path.join(project_path, f'products/full_disk/{filename}')
        corrected_map.save(img_savepath, overwrite=True)
        if not QUIET: print('FD saved:', filename)

def roi_correction(file):
    '''
    Applies correction to RoI images
    '''
    roi_map= Map(file)
    col, row= roi_map.meta['X1'], roi_map.meta['Y1']
    s_row, s_col= roi_map.meta['NAXIS1'], roi_map.meta['NAXIS2']
    roi_flat= flat_frame [row:row+s_row, col-20:col+s_col-20]
    corrected_roi_data= roi_map.data/roi_flat
    corrected_roi_data= np.nan_to_num(corrected_roi_data, nan=0.0)
    corrected_roi_data[corrected_roi_data> 6e4]=0
    corrected_roi_map= Map(corrected_roi_data, roi_map.meta)
    if SAVE:
        filename= roi_map.meta['F_NAME']
        img_savepath= os.path.join(project_path, f'products/roi/{filename}')
        corrected_roi_map.save(img_savepath, overwrite=True)
        if not QUIET: print('ROI saved:', filename)

if __name__=='__main__':
    SAVE= True # Toggle to save corrected image
    QUIET= False # Enable to run without output
    project_path= os.path.abspath('..')
    # Filepath for full disk images
    fd_files= sorted(glob.glob(os.path.join(project_path, 'data/raw/full_disk/*.fits'))) 
    # Filepath for ROI images
    roi_files= sorted(glob.glob(os.path.join(project_path, 'data/raw/roi/*.fits')))
    aligned_maps, unaligned_maps= alignmaps(fd_files) # Align maps
    flat_frame= makeflat(aligned_maps, unaligned_maps) # Make flat file. Global variable
    with multiprocessing.Pool() as pool:
        r1= pool.map_async(fd_correction, fd_files)
        if roi_files: # Check if roi_files is empty. Will run if not empty.
            r2= pool.map_async(roi_correction, roi_files)
            r1.get() # get() stops process if error is received. Use wait() for continued run, with errors.
            r2.get()
        else:
            r1.get()
