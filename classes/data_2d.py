from itertools import product

import numpy as np
from astropy.io import fits
import astropy.units as u
from myutils.logger import get_logger

from .data import Data

class Data2D(Data):
    """Defines a data in 2D.

    Can be used for loading images or visibility maps.
    
    Attributes:
        address: file name.
        data: the data.
        nhdu: HDU number to work with.
    """

    def __init__(self, address, nhdu=0):
        """Defines a new data object.

        Parameters:
            address (str): file name.
            nhdu (int, default=0): HDU number.
        """
        super(Data2D, self).__init__(address)
        self.nhdu = nhdu
        self.logger = get_logger(__name__)

        # Check dimensions
        if len(self.array.shape)>2:
            self.logger.warn('Reducing data dimensions: %r', self.array.shape)
            for i in range(len(self.array.shape)-2):
                self.array = self.array[0]
            self.logger.warn('New dimensions: %r', self.array.shape)

            # Update header
            self.header['NAXIS'] = 2
            keys = ['CDELT%i', 'CRVAL%i', 'CTYPE%i', 'CROTA%i', 'NAXIS%i']
            for i, key in product([3,4], keys):
                try:
                    self.logger.debug('Deleting %s', key % i)
                    del self.data[self.nhdu].header[key % i]
                except:
                    self.logger.debug('%s could not be deleted', key % i)
                    pass


    def load(self):
        """Load the data"""
        self.data = fits.open(self.address)

    def save(self, file_name=None):
        """Save the data.

        Parameters:
            file_name (default=None): new file name.
        """
        self.data.writeto(file_name or self.address , clobber=True)

    @property
    def array(self):
        return self.data[self.nhdu].data

    @array.setter
    def array(self, value):
        self.data[self.nhdu].data = value

    @property
    def header(self):
        return self.data[self.nhdu].header

    @header.setter
    def header(self, key, value):
        self.data[self.nhdu].header[key] = value

    @property
    def unit(self):
        try:
            return 1.*u.Unit(self.header['BUNIT'])
        except ValueError:
            return 1.*u.Unit(self.header['BUNIT'].lower())

    def max_pix(self):
        """Determine the position of the maximum.
        
        Zero based."""
        ymax, xmax = np.unravel_index(np.nanargmax(self.array),
            self.array.shape)
        return xmax, ymax
