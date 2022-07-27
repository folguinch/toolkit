import os

from astropy.stats import gaussian_fwhm_to_sigma
import astropy.units as u
import numpy as np

import toolkit.logger as logger

# Start settings
if os.path.isdir('logs'):
    LOG = logger.get_logger(__name__, filename='logs/myutils_mask.log')
else:
    LOG = logger.get_logger(__name__, filename='mytils_masks.log')

def circular_mask(array=None, shape=None, xy=None, ij=None, r=None, true_inside=True):
    # Check input
    if r is None:
        raise TypeError('Mask radius needed')
    if (xy is None and ij is None) or \
        (xy is not None and ij is not None):
        raise TypeError('Circular mask needs either xy or ij value')
    elif ij is not None:
        xy = ij[::-1]
    if (shape is None and array is None) or \
        (shape is not None and array is not None):
        raise TypeError('Circular mask needs either np.array or shape value')
    elif array is not None:
        if array.ndim != 2:
            raise ValueError('Array should be 2-D, shape = %r', array.ndim)
        else:
            shape = array.shape
    elif len(shape) != 2:
        raise ValueError('Shape should be 2-D, shape = %r', shape)

    # Distance
    yy, xx = np.indices(shape, dtype=float)
    dist = np.sqrt((yy - xy[1])**2 + (xx-xy[0])**2)

    # Mask
    mask = dist<=r
    if true_inside:
        LOG.info('Data in mask: %i', np.sum(mask))
        return mask
    else:
        LOG.info('Data in mask: %i', np.sum(~mask))
        return ~mask
        
def image_circular_mask(img, xy=None, ij=None, r=None, bmin=None, bmaj=None,
        true_inside=True):
    # Greet
    LOG.info('Calculating image circular mask')

    # Determine radius
    pixsize = np.sqrt(abs(img.header['CDELT1'] * img.header['CDELT2']))
    if r is not None:
        if hasattr(r, 'unit'):
            r = r.to(u.deg).value
            r = r / pixsize
        LOG.info('Input radius: %f pix', r)
    else:
        if bmin is None and bmaj is None:
            LOG.info('Using header bmin & bmaj')
            bmin = img.header['BMIN']
            bmaj = img.header['BMAJ']
        else:
            LOG.info('Using input bmin & bmaj')
        r = np.sqrt(bmin*bmaj) * gaussian_fwhm_to_sigma
        try:
            r = r.to(u.deg).value 
        except AttributeError:
            pass
        r = r / pixsize
        LOG.info('Mask radius: %f pix', r)

    # Use data shape
    shape = img.shape
    if len(shape)>2:
        shape = shape[-2:]

    return circular_mask(shape=shape, xy=xy, ij=ij, r=r,
            true_inside=true_inside)

