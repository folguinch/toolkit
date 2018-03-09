import numpy as np
from astropy.io import fits
import astropy.units as u
import scipy.ndimage as scimage
from scipy.stats import linregress

from .logger import get_logger
from .rotate import rotate, rotate_coord_clockwise

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

def add_stokes(img):
    """Add the Stokes axis.

    Parameters:
        fits (astropy.io.fits): fits file.
    """
    if img.data.ndim==4:
        return
    elif img.data.ndim==3:
        newfits = np.array([img.data])
        assert newfits.ndim==4
        newheader = img.header
        newheader['NAXIS'] = 4
        newheader['NAXIS4'] = 1
        newheader['CTYPE4'] = 'STOKES'
        newheader['CUNIT4'] = ''
        newheader['CRPIX4'] = 1
        newheader['CDELT4'] = 1
        newheader['CRVAL4'] = 1
        hdu = fits.PrimaryHDU(newfits)
        hdu.header = newheader
        return hdu
    else:
        newfits = np.array([img.data])
        assert newfits.ndim==3
        newheader = img.header
        newheader['NAXIS'] = 3
        newheader['NAXIS4'] = 1
        newheader['CTYPE3'] = 'STOKES'
        newheader['CUNIT3'] = ''
        newheader['CRPIX3'] = 1
        newheader['CDELT3'] = 1
        newheader['CRVAL3'] = 1
        hdu = fits.PrimaryHDU(newfits)
        hdu.header = newheader
        return hdu

def rotate_cube(img, angle, source, centre=(0,0), mode='bilinear', **kwargs):
    """Rotate image cube.

    Parameters
        angle: angle of rotation in degrees.
        centre: centre of rotation.
        mode: interpolation mode.
    """
    data = fits.open(img)[0]
    data = add_stokes(data)
    old_centre = data.data.shape[2]/2.-.5, data.data.shape[3]/2.-.5

    # Rotate the cube
    new_img = None
    naxis4 = int(data.header['NAXIS4'])
    for i in range(naxis4): 
        aux = None
        for j, slc in enumerate(data.data[i]):
            rotated = rotate(slc, angle, centre=centre, mode=mode)
            if aux is None:
                aux = np.array([rotated])
            else:
                aux = np.append(aux, [rotated], axis=0)
        if new_img is None:
            new_img = np.array([aux])
        else:
            new_img = np.append(new_img, [aux], axis=0)

    # Create a new fits file
    hdu = fits.PrimaryHDU(new_img)
    hdu.header = data.header
    hdu.header['COMMENT'] = 'Rotated %.3f deg, center %i, %i' % \
            ((angle,)+centre)
    
    # Redefine the reference pixel
    #if 'CRPIX1' in self.header:
    #    xref, yref = rotate_coord_clockwise(self['CRPIX1']-1,
    #                                        self['CRPIX2']-1, angle, 
    #                                        centre=old_centre[::-1])
    #    xref += new_img.shape[3]/2. - .5
    #    yref += new_img.shape[2]/2. - .5
    #    if kwargs.get('ref_pos'):
    #        ref_pos = kwargs['ref_pos']
    #    elif self.header.get('CRVAL1') and self.header.get('CRVAL2'):
    #        ref_pos = (self['CRVAL1'], self['CRVAL2'])
    #    else:
    #        ref_pos = None
    #    self.set_refpix(xref+1, yref+1, ref_pos)
    #if kwargs.get('centre_pos') is not None:
    #    self.set_centrepos(kwargs['centre_pos'])
    #elif kwargs.get('max_pos') is not None:
    #    self.set_maxpos(kwargs['max_pos'])
    #elif kwargs.get('set_pos') is not None:
    #    self.set_refpix(*kwargs['set_pos'])
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
