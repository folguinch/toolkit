from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from myutils.logger import get_logger

from .data import Data

class Data3D(Data):
    """ Creates a data cube structure.

    Attributes:
        address: file name.
        data: the data cube.
    """

    def __init__(self, address):
        """Defines a new data cube object.

        Parameters:
            address (str): file name.
        """
        super(Data3D, self).__init__(address)
        self.logger = get_logger(__name__)

    def load(self):
        """Load the data cube"""
        try:
            self.data = SpectralCube.read(self.address)
        except:
            self.logger.warn('Trying to fix the cube')
            cube = fits.open(self.address)[0]
            cube.header['CUNIT3'] = 'm/s'
            cube.header['CRVAL3'] = cube.header['CRVAL3'] * 1.E3
            cube.header['CDELT3'] = cube.header['CDELT3'] * 1.E3
            self.data = SpectralCube(cube.data, WCS(cube.header))
            self.logger.info('Cube fixed')

    def save(self, filename=None, overwrite=True):
        """Save the data cube"""
        self.data.write(filename or self.address, overwrite=overwrite)

    @property
    def wcs(self):
        """Return the WCS with the position data."""
        return self.data.wcs.sub(['longitude', 'latitude'])

    def get_coord(self, xpix, ypix, frame='fk5'):
        """Return the (ra, dec) coordinates of the input pixel location.

        Parameters:
            xpix (float): x-position of the coordinate.
            ypix (float): y-position of the coordinate.
            frame (str, default=fk5): sky frame projection.

        Returns:
            coord (astropy.SkyCoord): sky coordinate of the input location.

        Note:
            xpix and ypix are zero based.
        """
        ra, dec = self.wcs.all_pix2world([[xpix, ypix]], 0)[0]
        return SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame=frame)
