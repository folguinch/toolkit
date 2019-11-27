import numpy as np
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from ..logger import get_logger
from ..rotate import rotate
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

    def get_pixel(self, coord):
        crd = [[coord.ra.hour, coord.dec.degree]]
        return self.wcs.all_world2pix(crd, 0)[0]
        #return self.wcs.world_to_pixel(coord)

    def get_spectrum(self, coord=None, pixel=None):
        if coord is not None:
            x, y = self.get_pixel(coord)
        elif pixel is not None:
            x, y = pixel
        else:
            raise ValueError('Could not determine position')
        self.logger.info('Extracting spectrum at pixel: %i, %i', x, y)
        return self.data[:,y,x]

    def rotate(self, angle, source, centre=(0,0), mode='bilinear', **kwargs):
        """Rotate the cube.

        Parameters
            angle: angle of rotation in degrees.
            centre: centre of rotation.
            mode: interpolation mode.
        """
        #old_centre = self.centre

        # Rotate the cube
        aux = None
        for j, slc in enumerate(self.data.unmasked_data[:]):
            rotated = rotate(slc, angle, centre=centre, mode=mode)
            if aux is None:
                aux = np.array([rotated])
            else:
                aux = np.append(aux, [rotated], axis=0)

        # Create a new fits file
        hdu = fits.PrimaryHDU(aux)
        hdu.header = self.data.header
        hdu.header['COMMENT'] = 'Rotated %.3f deg, center %i, %i' % \
                ((angle,)+centre)
        self.logger.info('Image rotated %.1f deg, center %i, %i', angle, *centre)
        
        # Redefine the reference pixel
        hdu.header['CRPIX1'] = aux.shape[2]/2. + .5
        hdu.header['CRPIX2'] = aux.shape[1]/2. + .5
        hdu.header['CRVAL1'] = source.position.ra.to(u.deg).value
        hdu.header['CRVAL2'] = source.position.dec.to(u.deg).value

        # Redifine header rotation
        #if 'CROTA2' in hdu.header:
        #    #if hdu.header['CROTA2']==angle:
        #    self.logger.info('Deleting rotation keywords from header')
        try:
            del hdu.header['CROTA2']
            del hdu.header['CROTA1']
        except KeyError:
            pass
            #else:
            #    self.logger.info('Changing rotation keywords from header')
            #    hdu.header['CROTA2'] = hdu.header['CROTA2']-angle
        #else:
        #    self.logger.info('Setting rotation header keywords')
        #    hdu.header['CROTA1'] = 0
        #    hdu.header['CROTA2'] = angle

        return hdu
