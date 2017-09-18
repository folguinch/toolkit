import numpy as np
import astropy.units as u

"""Functions for working with numpy arrays"""

def load_struct_array(file_name, usecols=None):
    """Load a structured array table.

    Each file has a header. The first row has the name of each parameter and
    the second the units.

    Parameters:
        file_name (str): file to be loaded.
        usecols (iterable): columns to load.
    """

    with open(file_name, 'r') as input:
        # Read first two lines
        line1 = input.readline().strip(' #').split()
        line2 = input.readline().strip(' #').split()

        # Read dtype and data units
        dtype = []
        units = []
        for i in usecols:
            dtype += [(line1[i], float)]
            units += [u.Unit(line2[i])]

        # Read data
        data = np.loadtxt(input, usecols=usecols, dtype=dtype)

    return data, units
