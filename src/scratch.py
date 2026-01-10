import os
import glob
import matplotlib.pyplot as plt
from sunpy.map import Map
import astropy.units as u
import numpy as np
from scipy.ndimage import zoom
from sunkit_image.coalignment import apply_shifts
from astropy.convolution import convolve, Box2DKernel

def blur(data, kernel): #blurring function
    return(convolve(data, Box2DKernel(kernel), normalize_kernel=True))

def plot(array, VMX=None):
    plt.figure()
    plt.imshow(array, origin='lower', vmin=0, vmax=VMX)
    plt.show()

def get_submap(ref_img):
	"""
	Get a submap from the reference image for template matching.
	Parameters:
	- ref_img: Reference SunPy map image.
	Returns:
	- Submap of the reference image.
	"""
	center_coord = SkyCoord(0 * u.arcsec, 950* u.arcsec, frame=ref_img.coordinate_frame) #54,157
	width = 1100 * u.arcsec
	height =300 * u.arcsec   
	offset_frame = SkyOffsetFrame(origin=center_coord, rotation=0*u.deg)
	rectangle = SkyCoord(lon=[-1/2, 1/2] * width, lat=[-1/2, 1/2] * height, frame=offset_frame)
	ref_submap = ref_img.submap(rectangle) #bottom_left, top_right=top_right)
	return ref_submap

project_path= os.path.abspath('..')
# Filepath for full disk images
fd_files= sorted(glob.glob(os.path.join(project_path, 'data/raw/full_disk/*.fits'))) 
seq= Map(fd_files, sequence=True)
template_map=seq[0]
for m in seq:
    x_arry=[template_map.meta['CRPIX1']-m.meta['CRPIX1'] for i, m in enumerate(seq)]
    y_arry=[template_map.meta['CRPIX2']-m.meta['CRPIX2'] for i, m in enumerate(seq)]

#ref_submap = get_submap(template_map)
aligned_imgs= np.stack([m.data for m in aligned_maps], axis=0)
map_arr= np.stack([m.data for m in seq], axis=0)
raw_med= np.median(map_arr, axis=0)
med= np.median(aligned_imgs, axis=0)
med[med==0]=1

flat_frame= raw_med/med
flat_frame[flat_frame==0]=1
large_scale= blur(flat_frame, 10)
flat_frame=flat_frame/large_scale
corrected_img= template_map.data/flat_frame 
