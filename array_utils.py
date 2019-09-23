import numpy as np
import astropy.units as u

"""Functions for working with numpy arrays"""

def load_struct_array(file_name, usecols=None):
    """Load a structured array table.

    Each file has a header. The first row has the name of each parameter and
    the second the units. The name of each parameter must be different.

    Parameters:
        file_name (str): file to be loaded.
        usecols (iterable): columns to load.
    """
    with open(file_name, 'r') as input:
        # Read first two lines
        line1 = input.readline().strip(' #').split()
        line2 = input.readline().strip(' #').split()
        assert len(line1)==len(set(line1))

        # Read dtype and data units
        dtype = []
        units = {}
        for dty,unit in zip(line1,line2):
            dtype += [(dty, float)]
            units[dty] = 1.*u.Unit(unit)

        # Read data
        data = np.loadtxt(input, usecols=usecols, dtype=dtype)

    return data, units

def load_mixed_struct_array(file_name, usecols=None):
    """Load a structured array table of mixed types.

    Each file has a header. The first row has the name of each parameter and
    the second the units. The name of each parameter must be different. String
    values are marked with unit `-` in the second line.

    Parameters:
        file_name (str): file to be loaded.
        usecols (iterable): columns to load.
    """
    with open(file_name, 'r') as input:
        # Read first two lines
        line1 = input.readline().strip(' #').split()
        line2 = input.readline().strip(' #').split()
        assert len(line1)==len(set(line1))

        # Read dtype and data units
        names = []
        units = {}
        for i,(dty,unit) in enumerate(zip(line1,line2)):
            if usecols is not None and i not in usecols:
                continue
            names += [dty]
            if unit == '-':
                units[dty] = None
            else:
                units[dty] = 1.*u.Unit(unit)

        # Read data
        data = np.genfromtxt(input, usecols=usecols, dtype=None, names=names)

    return data, units

def save_struct_array(file_name, data, units, fmt='%10.4e\t'):
    r"""Save a structured array table.

    Save the data in a way it can be loaded by the load_struct_array function.

    Parameters:
        file_name (str): name of the file.
        data (nump.array): array to save.
        units (astropy.units dict): physical units of the data.
        fmt (str, default='10.4e\t'): string format for the data.
    """
    with open(file_name, 'w') as output:
        # Header
        line1 = '#' 
        line2 = '#' 
        for name in data.dtype.names:
            line1 += '{0}\t'.format(name)
            line2 += '{0}\t'.format(units[name])
        line1 = line1.strip() + '\n'
        line2 = line2.strip() + '\n'
        lines = line1 + line2

        # Data
        lines += '\n'.join((fmt*len(d)).strip() % tuple(d) for d in data)
        
        # Write
        output.write(lines)

        


