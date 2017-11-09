import numpy as np
import scipy.ndimage as scimage

"""Functions for working with fits files.
"""

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
    if axis==0:
        for col,ind in zip(img.T, mask.T):
            center += [np.sum(col[ind] * y)/np.sum(col[ind])]
    elif axis==1:
        for row,ind in zip(img, mask):
            center += [np.sum(row[ind]* x])/np.sum(row[ind])]
    else:
        center = scimage.measurements.center_of_mass(img)
    return np.array(center)

