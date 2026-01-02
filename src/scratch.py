import os
import glob
import matplotlib.pyplot as plt
from sunpy.map import Map
import astropy.units as u
import numpy as np
from sunkit_image.coalignment import apply_shifts
from astropy.convolution import convolve
from astropy.convolution import Box2DKernel

def blur(data, kernel): #blurring function
    return(convolve(data, Box2DKernel(kernel), normalize_kernel=True))

def plot(array):
    plt.figure()
    plt.imshow(array, origin='lower', vmin=0, vmax=1.5e4)
    plt.show()

project_path= os.path.abspath('..')
# Filepath for full disk images
fd_files= sorted(glob.glob(os.path.join(project_path, 'data/raw/full_disk/*.fits'))) 
seq= Map(fd_files, sequence=True)
template_map=seq[0]
for m in seq:
    x_arry=[template_map.meta['CRPIX1']-m.meta['CRPIX1'] for i, m in enumerate(seq)]
    y_arry=[template_map.meta['CRPIX2']-m.meta['CRPIX2'] for i, m in enumerate(seq)]

aligned_maps = apply_shifts(seq, yshift= y_arry * u.pixel, xshift=x_arry * u.pixel, clip=False)
aligned_imgs= np.stack([m.data for m in aligned_maps], axis=0)
med= np.median(aligned_imgs, axis=0)
med[med==0]=1
flat_frame= template_map.data/med
flat_frame[flat_frame==0]=1
large_scale= blur(flat_frame, 10)
flat_frame=flat_frame/large_scale
corrected_img= template_map.data/flat_frame 
