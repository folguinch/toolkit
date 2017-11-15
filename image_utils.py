import numpy as np
import astropy.units as u
import scipy.ndimage as scimage
from scipy.stats import linregress

"""Functions for working with fits files.
"""

def get_coord_axis(img, axis):
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
    pixsize = np.abs(img.header['CDELT%i' % axis]) * unit
    npix = img.header['NAXIS%i' % axis]

    return refval + (np.arange(npix, dtype=float) - (refpix-1))*pixsize

def get_coord_axes(img):
    """Get the coordinate axes from header.

    Parameters:
        img (astropy.hdu): input hdu.
    """
    x = get_coord_axis(img, 1)
    y = get_coord_axis(img, 2)

    return x, y

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

def lin_vel_gradient(img, sigma=0., nsigma=5):
    """Calculates the linear velocity gradient.

    Parameters:
        img (astopy.fits): input hdu.

    Returns:
        gradient (float): the velocity gradient.
        center (np.array): the center of mass.
    """
    assert hasattr(img, 'header')

    # Determine the axis for the center of mass
    if img.header['CTYPE1']=='OFFSET':
        axis = 1
        newx_units = u.Unit(img.header['CUNIT2'])
        newy_units = u.Unit(img.header['CUNIT1'])
    else:
        axis = 0
        newx_units = u.Unit(img.header['CUNIT1'])
        newy_units = u.Unit(img.header['CUNIT2'])

    # Calculate center of mass
    x, y = get_coord_axes(img)
    newx, center = center_of_mass(img.data, axis=axis, 
            mask=img.data>nsigma*sigma, x=x.value, y=y.value)
    print x.to(u.arcsec), y.to(u.km/u.s)

    # Fit line
    fit = linregress(center, newx)
    slope = fit.slope * newx_units / newy_units
    print slope.to(u.km/u.s/u.arcsec)

    return fit, center
