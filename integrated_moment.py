"""
Integrated moment map.
"""

from astronomy_image import AstroImage

class IntegratedMoment(AstroImage):
    """The moment 0 map of a radio observation."""
    def __init__(self, moment_0_name: str) -> None:
        """Initializing the moment map by reading in the file name."""
        super().__init__(moment_0_name)
        self.bmaj = self.header['BMAJ'] * 3600
        self.bmin = self.header['BMIN'] * 3600
        self.bpa = self.header['BPA']

    def create_beam_kernel(self, target_pix_scale: float):
        """Uses the beam to create a Gaussian 2d kernel"""
        return self.create_2d_gaussian_kernal((self.bmaj, self.bmin), target_pix_scale, self.bpa)

if __name__ == '__main__':
    infile = '/home/tlambert/Desktop/HZ7/Hz7_ISM/data/HZ7_mom_mask2.integrated.fits'
    hz7_mom0 = IntegratedMoment(infile)
