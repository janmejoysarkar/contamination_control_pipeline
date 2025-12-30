#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sun Oct 19 08:02:45 PM CEST 2025
@author: janmejoy
@hostname: machine

DESCRIPTION
"""

import numpy as np
from astropy.io import fits
from astropy.convolution import Box2DKernel
import os
import glob
import matplotlib.pyplot as plt
from skimage.filters import median

project_path= os.path.abspath('..')
data_path= os.path.join(project_path,'data/raw/SUT_T25_1555_001464_Lev1.0_2025-10-19T19.45.33.559_0971NB03.fits')
VMX= 1.3e4
SAVE=False
savepath= os.path.join(project_path,'reports/savefig.png')
row, col, size= 1700, 2550, 200

#kernel= np.array([[1]*10]) # 10 px horz
#kernel= np.array([[1*10]]) # 10 px vert
kernel= Box2DKernel(10) #10 px box

with fits.open(data_path) as hdu:
    data= hdu[0].data
    head= hdu[0].header

med= median(data, footprint=kernel)

plt.figure()
plt.subplot(2,2,1)
plt.imshow(data, vmin=0, vmax=VMX, origin='lower')
plt.title('raw')
plt.subplot(2,2,2)
plt.imshow(med, vmin=0, vmax=VMX, origin='lower')
plt.title('filtered')
plt.subplot(2,2,3)
plt.imshow(data[row:row+size, col:col+size], vmin=0, vmax=VMX, origin='lower')
plt.subplot(2,2,4)
plt.imshow(med[row:row+size, col:col+size], vmin=0, vmax=VMX, origin='lower')
if SAVE:
    plt.savefig(savepath, dpi=300)
    plt.close()
else:
    plt.show()
