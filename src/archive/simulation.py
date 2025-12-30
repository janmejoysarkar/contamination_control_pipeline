#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2025-07-22 17:21:41
@author: janmejoyarch
@hostname: suitpoc1

DESCRIPTION

- Image generator to simulate contamination spots
- Prepared to test KLL flat fielding.
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os


def make_circ(r,c,rad,val):
    #returns a mask
    X,Y= np.meshgrid(x,x)
    circ= np.sqrt((X-c)**2+(Y-r)**2)
    circ= circ < rad
    circ=circ*val
    return circ

SIZE=4096 # size of base square image
RSUN=1400 # px
SUN_COUNT=1e4
BIAS=500
DITH=30 # +- pixels
N=20
SAVE= True

x= np.arange(SIZE)

project_path= os.path.abspath('..')

### Dust Model ###
dust_model= np.ones((SIZE,SIZE))
REPS= 50
d_xs= np.random.randint(SIZE//2-RSUN, SIZE//2+RSUN, size=REPS)
d_ys=np.random.randint(SIZE//2-RSUN, SIZE//2+RSUN, size=REPS)
vals= np.random.uniform(0.4, 0.8, size=REPS)
radii=np.random.randint(7,30,size=REPS)

for d_x, d_y, val, rad in zip(d_xs, d_ys, vals, radii):
    dust= make_circ(d_x, d_y, rad, val)
    dust[dust==0]=1
    dust_model= dust_model*dust

### Sun Images ###
row_dithers= np.random.randint(-DITH, DITH, size=N)
col_dithers= np.random.randint(-DITH, DITH, size=N)
imgs= []
imgname=1
for row_dither, col_dither in zip(row_dithers, col_dithers):
    sun=make_circ(SIZE//2+col_dither, SIZE//2+row_dither, RSUN, SUN_COUNT)
    img= sun*dust_model
    img= img+BIAS # add bias
    imgs.append(img)
    if SAVE:
        fits.writeto(os.path.join(project_path, f'data/interim/img_{imgname:02d}.fits'), img, overwrite=True)
        imgname+=1

if SAVE:
    fits.writeto(os.path.join(project_path, f'data/interim/dust_model.fits'), dust_model, overwrite=True)
    dither_array= np.column_stack((row_dithers, col_dithers))
    np.savetxt(os.path.join(project_path, f'data/interim/dither.txt'), dither_array, header='Row_Dith, Col_Dith')
