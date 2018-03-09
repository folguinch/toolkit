import scipy.ndimage.interpolation as interp
import scipy.interpolate as interp2
import numpy as np

from .interpolation import bilinear_interp

def new_round(x):
    fl = np.floor(x)
    mask = x-fl >= 0.5
    try:
        fl[mask] = fl[mask]+1
    except IndexError:
        if mask:
            fl += 1
    
    return fl

def rotate_coord_counter(x, y, angle, centre=[0, 0]):
    """
    Obtain the new coordinates of the points (x,y) after a counter-clockwise 
    rotation by a given angle around a given rotation centre.

    Parameters:
    -----------
    x, y: float or ndarray
        x-axis and y-axis coordinates
    angle: float
        rotation angle in degrees
    centre: array
        position of the rotation centre
    """
    theta = np.radians(angle)
    newx = (x-centre[0])*np.cos(theta) - (y-centre[1])*np.sin(theta)
    newy = (x-centre[0])*np.sin(theta) + (y-centre[1])*np.cos(theta)
    return newx, newy

def rotate_coord_clockwise(x, y, angle, centre=[0, 0]):
    """
    Obtain the new coordinates of the points (x,y) after a clockwise
    rotation by a given angle around a given rotation centre.
    
    Parameters:    
    -----------
    x, y: float or ndarray
        x-axis and y-axis coordinates
    angle: float
        rotation angle in degrees
    centre: array
        position of the rotation centre
    """
    theta = np.radians(angle)
    newx = (x-centre[0])*np.cos(theta) + (y-centre[1])*np.sin(theta)
    newy = -(x-centre[0])*np.sin(theta) + (y-centre[1])*np.cos(theta)
    return newx, newy

def rotate(im, ang, centre=(0,0), mode='bilinear'):

    n90 = int(ang//90)
    if n90 != 0:
        img = np.rot90(im, n90)
        angle = ang - 90.*n90
        if angle%90 == 0:
            return img
    else:
        angle = ang
        img = im

    if mode=='bilinear':
        return rotate_bilinear(img, angle, centre=centre)
    elif mode=='nearest':
        return rotate_nearest(img, angle, centre=centre)

def rotate_bilinear(img, angle, centre=(0,0)):
    X_cen, Y_cen = rotate_coord_clockwise(centre[0], centre[1], angle)

    xsize = np.abs(img.shape[1]*np.cos(np.radians(angle))) + \
            np.abs(img.shape[1]*np.cos(np.radians(angle+90)))
    ysize = np.abs(img.shape[0]*np.sin(np.radians(angle))) + \
            np.abs(img.shape[0]*np.sin(np.radians(angle+90)))
    xsize = int(round(ysize))
    ysize = int(round(xsize))
    img_new = np.ones((ysize, xsize))*np.nan
    Y_new, X_new = np.indices(img_new.shape)

    X_new_rot, Y_new_rot = rotate_coord_counter(X_new, Y_new, angle,
                                                centre=[X_cen, Y_cen])
    h = img.shape[0] * np.sin(np.radians(angle))
    dx = h * np.cos(np.radians(90.-angle)) + centre[0]
    dy = h * np.sin(np.radians(90.-angle)) - centre[1]
    X_new_rot += dx
    Y_new_rot -= dy
    mask_x = (X_new_rot >= 0) & (X_new_rot <= img.shape[1])
    mask_y = (Y_new_rot >= 0) & (Y_new_rot <= img.shape[0])
    mask = mask_x & mask_y

    img_new[mask] = bilinear_interp(img, X_new_rot[mask], Y_new_rot[mask])

    return img_new


def rotate_nearest(img, ang, centre=(0,0)):

    Y, X = np.indices(img.shape)
    
    X_new = X*np.cos(np.radians(ang)) + Y*np.sin(np.radians(ang))
    Y_new = -X*np.sin(np.radians(ang)) + Y*np.cos(np.radians(ang))
    
    X_new = new_round(X_new).astype(int) #+ np.min(X_new)
    Y_new = new_round(Y_new).astype(int) #+ np.min(Y_new)

    xsize = np.abs(img.shape[1]*np.cos(np.radians(ang))) + \
        np.abs(img.shape[1]*np.cos(np.radians(ang+90)))
    ysize = np.abs(img.shape[0]*np.sin(np.radians(ang))) + \
        np.abs(img.shape[0]*np.sin(np.radians(ang+90)))
    xsize = int(round(ysize))
    ysize = int(round(xsize))
    img_new = np.zeros((ysize, xsize))
    img_new[np.ravel(Y_new),np.ravel(X_new)] = img[np.ravel(Y),np.ravel(X)]

    # Re-calculate zeros interpolating
    zeros = (img_new==0) & (np.roll(img_new,1,axis=1)!=0) & \
            (np.roll(img_new,-1,axis=1)!=0) & (np.roll(img_new,1,axis=0)!=0) & \
            (np.roll(img_new,-1,axis=0)!=0)
    Y, X = np.indices(img_new.shape)
    m = (np.roll(img_new, 1, axis=1)-np.roll(img_new, -1, axis=1)) / \
            (np.roll(X, 1, axis=1)-np.roll(X, -1, axis=1))
    n = np.roll(img_new, 1, axis=1) - np.roll(X, 1, axis=1)*m
    img1 = X*m + n
    m = (np.roll(img_new, 1, axis=0)-np.roll(img_new, -1, axis=0)) / \
            (np.roll(Y, 1, axis=0)-np.roll(Y, -1, axis=0))
    n = np.roll(img_new, 1, axis=0) - np.roll(Y, 1, axis=0)*m
    img2 = Y*m + n
    img3 = (img1+img2)/2.
    img_new[zeros] = img3[zeros]

    # Check new size
    if img.shape[0] == img.shape[1] and img_new.shape[0]!=img_new.shape[1]:
        raise Exception('Wrong image size')

    return np.roll(np.roll(img_new, np.abs(np.min(Y_new)), axis=0),
                   np.abs(np.min(X_new)), axis=1)
