import numpy as np
import astropy.units as u
import scipy.ndimage as scimage
from scipy.stats import linregress

from .logger import get_logger

"""Functions for working with fits files.
"""

def get_coord_axis(img, axis, pixsize=None):
    """Get the coordinate axis.
    
    Parameters:
        img (astropy.hdu): hdu.
        axis (int): axis to get (1 for x or 2 for y).
    """
    assert axis in [1, 2]
    assert hasattr(img, 'header')

    # Parameters from header
    unit = u.Unit(img.header['CUNIT%i' % axis])
    refval = img.header['CRVAL%i' % axis] * unit
    refpix = img.header['CRPIX%i' % axis]
    pixsize = pixsize or np.abs(img.header['CDELT%i' % axis]) * unit
    npix = img.header['NAXIS%i' % axis]

    # Check that the units are the same
    refval = refval.to(pixsize.unit)

    return refval + (np.arange(npix, dtype=float) - (refpix-1))*pixsize

def get_coord_axes(img, pixsizes=(None,None)):
    """Get the coordinate axes from header.

    Parameters:
        img (astropy.hdu): input hdu.
    """
    return [get_coord_axis(img, i, pixsize=pixsizes[i-1]) for i in [1,2]]


def center_of_mass(img, axis=None, mask=None, x=None, y=None):
    """Calculates the center of mass of img.

    The center of mass is calculated as:
           sum(img[i]*r[i])
      R =  ----------------
              sum(img[i])
    If axis is specified then r[i] are replaced by the values along that axis,
    e.g. if axis=0 then r[i]=y[i]. Note that axis uses the numpy standard, i.e.
    y-axis is 0 and x-axis is 1.

    Parameters:
        img (np.array): the image data
        axis (int, default=None): axis along which the weighted average is
            taken.
        mask (bool, default=None): mask for img. Only points in *img[mask]* are
            used.
        x (np.array, default=None): values of the x-axis.
        y (np.array, default=None): values of the y-axis.

    Returns:
        center (np.array): the value of the center of mass. A coordinate in
            terms of x and y if axis is not specified or an array of points if
            axis is specified.
    """
    # Validate input
    assert len(img.shape)==2
    assert axis in [0, 1, None]
    assert mask is None or mask.shape==img.shape
    assert x is None or len(x)==img.shape[1]
    assert y is None or len(y)==img.shape[0]

    # Define axes
    if x is None:
        x = np.arange(img.shape[1], dtype=float)
    if y is None:
        y = np.arange(img.shape[0], dtype=float)

    center = []
    newx = []
    if axis==0:
        for i,(col,ind) in enumerate(zip(img.T, mask.T)):
            if np.sum(ind)<=1:
                continue
            center += [np.sum(col[ind] * y[ind])/np.sum(col[ind])]
            newx += [x[i]]
    elif axis==1:
        for i,(row,ind) in enumerate(zip(img, mask)):
            if np.sum(ind)<=1:
                continue
            center += [np.sum(row[ind] * x[ind])/np.sum(row[ind])]
            newx += [y[i]]
    else:
        center = scimage.measurements.center_of_mass(img)
        return center

    return np.array(newx), np.array(center)

def lin_vel_gradient(img, sigma=0., nsigma=5, pixsizes=(None,None), 
        filterx=None, logger=get_logger(__name__)):
    """Calculates the linear velocity gradient.

    Parameters:
        img (astopy.fits): input hdu.

    Returns:
        gradient (float): the velocity gradient.
        center (np.array): the center of mass.
    """
    assert hasattr(img, 'header')
    assert len(pixsizes)==2
    assert filterx is None or len(filterx) in [1,2]

    # Determine the axis for the center of mass
    if img.header['CTYPE1']=='OFFSET':
        axis = 1
        if None not in pixsizes:
            newx_units = pixsizes[1].unit
            newy_units = pixsizes[0].unit
        else:
            newx_units = u.Unit(img.header['CUNIT2'])
            newy_units = u.Unit(img.header['CUNIT1'])
    else:
        axis = 0
        if None not in pixsizes:
            newx_units = pixsizes[0].unit
            newy_units = pixsizes[1].unit
        else:
            newx_units = u.Unit(img.header['CUNIT1'])
            newy_units = u.Unit(img.header['CUNIT2'])

    # Calculate center of mass
    x, y = get_coord_axes(img, pixsizes=pixsizes)
    mask = np.ones(img.data.shape, dtype=bool)
    if filterx is not None:
        Y,X = np.indices(img.data.shape)
        mask = X>filterx[0]
        if len(filterx)==2:
            mask = mask & (X<filterx[1])
    mask = mask & (img.data>nsigma*sigma)

    newx, center = center_of_mass(img.data, axis=axis, mask=mask, 
            x=x.value, y=y.value)

    # Fit line
    fit = linregress(center, newx)
    slope = fit.slope * newx_units / newy_units
    stderr = fit.stderr * newx_units / newy_units
    intercept = fit.intercept * newy_units
    logger.info('Linear velocity gradient fit:')
    logger.info('\tSlope = %s', slope.to(u.km/u.s/u.arcsec))
    logger.info('\tError = %s', stderr.to(u.km/u.s/u.arcsec))

    return slope, stderr, intercept, newx*newx_units, center*newy_units
