"""
Helper functions.
"""

import numpy as np
from astropy.io.fits.header import Header
from astropy.io.fits.hdu.hdulist import HDUList
from astropy.io.fits.hdu.image import ImageHDU

def determine_hdu_index(hdul: HDUList) -> int:
    """Determines the hdu index which has the correct data and header."""
    image_hdu = False
    for i, hdu in enumerate(hdul):
        if isinstance(hdu, ImageHDU):
            image_hdu = True
            index = i

    if not image_hdu:
        for i, hdu in enumerate(hdul):
            try:
                _  = hdul[0].data
                index = i
                break
            except AttributeError:
                pass

    return index

def read_pixscale_from_header(header: Header) -> float:
    """Determines the pixel scale from the given header object."""
    try:
        pix_scale = np.abs(header['CDELT2'])
        if pix_scale == 1:
            pix_scale = header['CD2_2']
    except KeyError:
        pix_scale = header['CD2_2']
    return pix_scale
