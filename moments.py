"""
Radio data class
"""

import numpy as np
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
import helper_funcs as hf

class IntegratedMoment:
    """The moment 0 map of a radio observation."""
    def __init__(self, cube_name: str) -> None:
        """Initializing the moment map by reading in the file name."""
        self.hdul = fits.open(cube_name)
        self.bmaj = self.hdul[0].header['BMAJ']
        self.bmin = self.hdul[0].header['BMIN']
        self.bpa = self.hdul[0].header['BPA']
    

if __name__ == '__main__':
    infile = '/home/tlambert/Desktop/HZ7/Hz7_ISM/data/HZ7_mom_mask2.integrated.fits'
    hz7_mom0 = IntegratedMoment(infile)
