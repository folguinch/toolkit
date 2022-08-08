"""Maths tools."""
from typing import Tuple, Optional

import astropy.stats as apystats
import numpy as np

def rms(x: np.array):
    """Root mean square (rms).

    The rms of `x` is defined as:
    ```
    rms = sqrt( sum(x_i**2) / n )
    ```
    where `n` is the number of points in `x`.

    Args:
      x: data.

    Returns:
      The root mean square of `x`.
    """
    return np.sqrt(np.sum(x**2)/x.size)

def quick_rms(data: np.array):
    """A quick estimation of the rms of the data."""
    return apystats.mad_std(data, ignore_nan=True)

def distance_array(shape: Tuple, position: Tuple,
                   mask: Optional[np.array] = None):
    """Returns an array with the distance of each point to the position.

    The `position` corresponds to the `(x, y)` position, i.e. in Cartesian
    coordinates.

    Args:
      shape: shape of the distance matrix.
      position: central point.
      mask: optional; masked distance value (`True` if masked).
    """
    ymesh, xmesh = np.indices(shape)
    dist = np.sqrt((xmesh - position[0])**2 + (ymesh - position[1])**2)

    if mask is not None:
        dist = np.ma.array(dist, mask=mask)

    return dist
