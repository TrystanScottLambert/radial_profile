"""
Helper functions.
"""

import numpy as np
from astropy.io.fits.hdu.table import BinTableHDU
from astropy.io.fits.hdu.hdulist import HDUList


def determine_average_beam_properties(beam_table: BinTableHDU) -> list[float]:
    """
    Reads in a beam table and determines the
    average BMAJ, BMIN, and BPA.
    """

    avg_bmaj = np.mean(beam_table.data.field('BMAJ'))
    avg_bmin = np.mean(beam_table.data.field('BMIN'))
    avg_bpa = np.mean(beam_table.data.field('BPA'))

    return avg_bmaj, avg_bmin, avg_bpa

def find_index_of_beam_table(hdul: HDUList) -> int:
    """
    Searches a hdu list and determines
    the position of the beam table if it exists.
    """
    for i, hdu_item in enumerate(hdul):
        if isinstance(hdu_item, BinTableHDU):
            if 'BMAJ' in str(hdu_item.header):
                return i
    print('No beam table found.')
