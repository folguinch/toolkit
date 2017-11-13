import numpy as np

def rms(x):
    """Root mean square (rms).

    The rms of *x* is defined as:
    rms = sqrt( sum(x_i**2) / n )
    where n is the number of points in x.

    Parameters:
        x (np.array): data.

    Return:
        rms (float): the root mean square of *x*.
    """
    assert hasattr(x, 'size')
    return np.sqrt(np.sum(x**2)/data.size)
