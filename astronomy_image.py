"""
Base class for 2D Astronomy images
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
import helper_funcs as hf

class AstroImage:
    """Base class for 2-dimensional astronomy images."""
    def __init__(self, file_name: str) -> None:
        """ Only requires the file name."""
        self.hdul = fits.open(file_name)
        self.hdu_index = hf.determine_hdu_index(self.hdul)
        self.header = self.hdul[self.hdu_index].header
        self.data = self.hdul[self.hdu_index].data
        self.wcs = WCS(self.header)
        self.pix_scale = hf.read_pixscale_from_header(self.header)
        self.rms = sigma_clipped_stats(self.data, mask = np.isnan(self.data))[2]

    def convolve(self, kernel: Gaussian2DKernel) -> np.ndarray:
        """Convolves the current data with a given 2d gaussin kernel."""
        convolved_data = convolve_fft(self.data, kernel)
        return convolved_data

    def create_2d_gaussian_kernal(
            self, psf: tuple[float, float],
            pix_scale:float, position_angle: float) -> Gaussian2DKernel:
        """
        Creates a 2d gaussian kernal which can be used to convolve other images.
        The given psf would be equivalent to the beam in a radio image and would then be (bmaj, bmin).
        For a HST/JWST like image the psf would be the 2D gaussian psf, (FWHM_maj, FWHM_min) but 
        this can be the (FWHM) where the psf if just a simple 1D gaussin.

        pix_scale is the pixel scale of the target image, not of the current image.

        position_angle is the position angle in degrees.
        """

        gaussian_sigma_2d = [(axes/gaussian_fwhm_to_sigma)/pix_scale for axes in psf]
        kernel = Gaussian2DKernel(x_stddev = gaussian_sigma_2d[0], y_stddev = gaussian_sigma_2d[1],
                                        theta = position_angle * (np.pi/180),
                                        x_size = self.data.shape[0], y_size = self.data.shape[1])
        return kernel

if __name__=='__main__':
    infile_hst = '/home/tlambert/Desktop/HZ7/Hz7_ISM/HZ7_Analysis/data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits'
    infile_alma = '/home/tlambert/Desktop/HZ7/Hz7_ISM/data/HZ7_mom_mask2.integrated.fits'

    hst = AstroImage(infile_hst)
    alma = AstroImage(infile_alma)