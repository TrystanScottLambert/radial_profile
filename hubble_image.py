"""
HST image
"""

from astronomy_image import AstroImage

class HubbleImage(AstroImage):
    """
    Main class for hubble images.
    """

if __name__ == '__main__':
    infile = '/home/tlambert/Desktop/HZ7/Hz7_ISM/HZ7_Analysis/data/ASCImages/gaiacorrected/hst_13641_07_wfc3_ir_f105w_sci_gaia_corrected.fits'
    hst_105 = HubbleImage(infile)
