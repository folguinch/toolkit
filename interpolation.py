import numpy as np

def bilinear_interp(img, x, y):
    """
    Notes:
        From Wikipedia 'Bilinear interpolation' and
        http://stackoverflow.com/questions/12729228/simple-efficient-bilinear-\
            interpolation-of-images-in-numpy-and-python
    """
    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, img.shape[1]-1)
    x1 = np.clip(x1, 0, img.shape[1]-1)
    y0 = np.clip(y0, 0, img.shape[0]-1)
    y1 = np.clip(y1, 0, img.shape[0]-1)

    I00 = img[y0, x0]
    I01 = img[y1, x0]
    I10 = img[y0, x1]
    I11 = img[y1, x1]

    f00 = (x1-x) * (y1-y)
    f10 = (x-x0) * (y1-y)
    f01 = (x1-x) * (y-y0)
    f11 = (x-x0) * (y-y0)

    return f00*I00 + f01*I01 + f10*I10 + f11*I11
