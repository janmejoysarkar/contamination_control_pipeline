#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2025-07-08 22:44:22
@author: janmejoyarch
@hostname: suitpoc1

DESCRIPTION
"""
from astropy.io import fits
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from scipy.ndimage import shift, gaussian_filter
from skimage.registration import phase_cross_correlation
from scipy.signal import correlate2d

proj_path= os.path.abspath('..')
datapath_list= glob.glob(os.path.join(proj_path,'data/raw/*'))

def read_data(datapath_list):
    crpix1_ls=[]
    crpix2_ls=[]
    data_ls=[]
    for datapath in datapath_list:
        with fits.open(datapath) as hdu:
            data=hdu[0].data
            hdr=hdu[0].header
            crpix1, crpix2= hdr['CRPIX1'], hdr['CRPIX2']
            data_ls.append(data)
            crpix1_ls.append(int(round(crpix1,0)))
            crpix2_ls.append(int(round(crpix2,0)))
    return(data_ls)

data_list= read_data(datapath_list)
image1, image2= data_list[1], data_list[9]
image1= gaussian_filter(image1/np.max(image1), sigma=5)
image2= gaussian_filter(image2/np.max(image1), sigma=5)
shift, error, phasediff= phase_cross_correlation(image1, image2, upsample_factor=10)
#image2= shift(image2, shift=(-10, -10), order=1)
#plt.figure()
#plt.imshow(image1-image2, origin='lower')
#plt.show()
