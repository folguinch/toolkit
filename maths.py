"""Maths tools."""
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
