#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2025-07-23 15:44:28
@author: janmejoyarch
@hostname: suitpoc1

DESCRIPTION
"""
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import shift
import numpy as np
import glob
import os

def save_jpg(imgs):
    i=0
    for img in imgs:
        plt.figure(figsize=(6, 6))
        plt.axis('off')
        plt.imshow(img, cmap='gray', origin='lower')
        plt.savefig(f'/home/janmejoyarch/Desktop/{i:02d}.jpg', dpi=300, bbox_inches='tight', pad_inches=0)
        print(f'{i:02d}')
        i+=1
        plt.close()

project_path= os.path.abspath('..')
img_files= sorted(glob.glob(os.path.join(project_path, 'data/interim/img*.fits')))
dither_file=os.path.join(project_path, 'data/interim/dither.txt')

imgs=[]
for img_file in img_files:
    with fits.open(img_file) as hdu:
        data= hdu[0].data
        imgs.append(data)

dither_data= np.loadtxt(dither_file)
dither_data= dither_data-dither_data[0]

## Image alignment wrt first image ##
aligned_imgs=[]
for i in range(len(imgs)):
    aligned_img= shift(imgs[i],(-dither_data[i][1], -dither_data[i][0]), order=3, mode='constant', cval=500)
    aligned_imgs.append(aligned_img)

## KLL ##
ff_ls=[]
s_ls=[]

num_iterations=10
# Initialize flat field
ff = np.ones_like(aligned_imgs[0])

# Outer loop: iterate KLL algorithm
for k in range(num_iterations):
    s_stack = []

    # Step 1: Estimate solar images from current flat field
    for img in aligned_imgs:
        s_est = img / ff
        s_stack.append(s_est)

    S = np.mean(s_stack, axis=0)  # Average solar scene

    ff_stack = []

    # Step 2: Estimate flat field from solar scene
    for img in aligned_imgs:
        ff_est = img / S
        ff_stack.append(ff_est)

    ff = np.mean(ff_stack, axis=0)  # Average flat field

    # Optional: normalize, smooth, display
    ff /= np.mean(ff)
    plt.imshow(S, cmap='gray'); plt.title(f'Iteration {k}'); plt.show()


'''
for i in np.arange(0,10):
    s= aligned_imgs[i]/ff#update the s
    ff= aligned_imgs[i]/s
    plt.imshow(s); plt.show()
s0= aligned_imgs[0]/ff0
ff1= aligned_imgs[1]/s0
s1= aligned_imgs[1]/ff1
ff2= aligned_imgs[2]/s1
s2= aligned_imgs[2]/ff2
ff3= aligned_imgs[3]/s2
s3= aligned_imgs[3]/ff3
ff4= aligned_imgs[4]/s3
s4= aligned_imgs[4]/ff4
for i in np.arange(1,10):
    s= shift(s, (-dither_data[i][1], -dither_data[i][0]), order=3, mode='constant', cval=500)
    ff= imgs[i]/s
    s= imgs[i]/ff
    print('iter:', i)
    ff_ls.append(ff)
    s_ls.append(s)
    #plt.imshow(s-imgs[i]); plt.show()
s= imgs[0]/ff
s= shift(s, (dither_data[1][1], dither_data[1][0]), order=3, mode='constant', cval=500)
ff= imgs[1]/s
ff= shift(ff, (-dither_data[2][1], dither_data[2][0]), order=3, mode='constant', cval=500)

s= imgs[2]/ff
s= shift(s, (dither_data[3][1], dither_data[3][0]), order=3, mode='constant', cval=500)
ff= imgs[3]/s
ff= shift(ff, (-dither_data[4][1], dither_data[4][0]), order=3, mode='constant', cval=500)
s= imgs[4]/ff
'''
