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
        