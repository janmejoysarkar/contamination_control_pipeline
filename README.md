
![Logo](https://suit.iucaa.in/sites/default/files/top_banner_compressed_2_1.png)


# SUIT Contamination Correction 🧹 ☀️ 🛰️ 

The SUIT instrument on Aditya-L1 maintains its image sensor at a very low temperature (-55 degC). This causes some volatiles to condense on the CCD surface. Traces of these contaminants are seen in the images recorded by SUIT.
This collection of modules is necessary to remove the contaminants from the SUIT images.

## Principle
The Aditya-L1 satellite has a jitter about its three axes. Jitter about the Pitch and Roll axes causes the image to shift randomly on the SUIT CCD. We exploit this as a dither by computing the exact shifts in the image, accurate to the pixel.

When the solar features are aligned, one particular contaminant pattern covers different regions of the Sun. This is exploited to generate a flat-field calibration file to remove the contaminants from each un aligned solar image.

## Modules 
The modules are designed to be user friendly- for implementation by the end user.

### `continuum_correction.py`: 
- Applies contaminant correction on full disk continuum images.

### `2k_fulldisk_correction.py`
- Used to correct 2k line channel full disk images.
- Uses sun-center information for co-alignment.
- Ensure the first 10 images of the input list are taken within a short time interval.

### `line_fulldisk_correction.py`
- Uses north limb of the sun for co-alignment.
- Useful for NB03 and NB04 channel images.

## `scratch.py`
- Scratchpad for testing code sections.

## Authors

- [@janmejoysarkar](https://github.com/janmejoysarkar)

## Acknowledgements

 - [ISRO, Aditya-L1](https://www.isro.gov.in/Aditya_L1.html)
## Screenshots

![Correction of contaminants on SUIT NB05 filter](README_files/Figure_1.png)


## Usage/Examples
SUIT image files are to be symlinked or placed at the specific folders based on the file type.

Folder structure for data

    ./data/raw/roi
    ./data/raw/full_disk

Folder structure for products

    ./products/roi
    ./products/full_disk


