import os, argparse
from configparser import ExtendedInterpolation

import astropy.coordinates as apycoord
from astropy.io import fits
import astropy.units as u
from astropy import wcs
import numpy as np

from ..myconfigparser import myConfigParser
from ..classes.dust import Dust
from ..logger import get_logger

def validate_files(filenames, check_is_file=True):
    try:
        validated = os.path.expanduser(filenames)
        if check_is_file and not os.path.isfile(validated):
            raise IOError('%s does not exist' % validated)
    except AttributeError:
        validated = []
        for fname in filenames:
            validated += [os.path.expanduser(fname)]
            if check_is_file and not os.path.isfile(validated[-1]):
                raise IOError('%s does not exist' % (validated[-1]))

    return validated

##### Loaders #####

class LoadConfig(argparse.Action):
    """Action class for loading a configuration file in argparse"""

    def __call__(self, parser, namespace, values, option_string=None):
        assert os.path.isfile(values)
        config = myConfigParser(interpolation=ExtendedInterpolation())
        config.read(values)
        setattr(namespace, self.dest, config)

class LoadArray(argparse.Action):
    """Action class for loading a np.array from command line"""

    def __call__(self, parser, namespace, values, option_string=None):
        array = np.array(values, dtype=float)
        setattr(namespace, self.dest, array)

class ArrayFromRange(argparse.Action):
    """Action class for creating a np.array with linspace from command line"""

    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        if nargs not in range(2, 5):
            raise ValueError('only 2,3 or 4 nargs allowed')
        super(ArrayFromRange, self).__init__(option_strings, dest,
                nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        start, stop = map(float, values[:2])
        if len(values)==4:
            n = int(values[2])
            if values[-1]=='linear':
                fn = np.linspace
            elif values[-1]=='log':
                fn = np.logspace
                start = np.log10(start)
                stop = np.log10(stop)
            else:
                raise NotImplementedError('%s not implemented' % values[-1])
        elif len(values)==3:
            fn = np.linspace
            n = int(values[2])
        else:
            fn = np.linspace
        value = fn(start, stop, n)
        setattr(namespace, self.dest, value)

class LoadStructArray(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(LoadStructArray, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from .array_utils import load_struct_array
        array = load_struct_array(validate_files(values))
        setattr(namespace, self.dest, array)

class LoadMixedStructArray(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(LoadMixedStructArray, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from .array_utils import load_mixed_struct_array
        array = load_mixed_struct_array(validate_files(values))
        setattr(namespace, self.dest, array)

class LoadTXTArray(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        array = np.loadtxt(values, dtype=float)
        setattr(namespace, self.dest, array)

class LoadFITS(argparse.Action):
    """Action for loading a FITS file with astropy"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values)
        try:
            vals = fits.open(''+values)[0]
        except TypeError:
            vals = []
            for val in values:
                vals += [fits.open(val)[0]]
        setattr(namespace, self.dest, vals)

class LoadDust(argparse.Action):
    """Action for loading dust files"""

    def __call__(self, parser, namespace, values, option_string=None):
        dust = Dust(values)
        setattr(namespace, self.dest, dust)

class LoadTable(argparse.Action):
    """Action for loading astropy Tables"""

    def __call__(self, parser, namespace, values, option_string=None):
        from .tables import Table

        try:
            tabname = ''+values
            #table_id = os.path.splitext(os.path.basename(tabname))[0]
            table = Table(tabname)
        except TypeError:
            if len(values)==2:
                table = Table(values[0], table_id=values[1])
            elif len(values)==1:
                tabname = values[0]
                #table_id = os.path.splitext(os.path.basename(tabname))[0]
                table = Table(tabname)
            else:
                raise ValueError('Number of values not allowed.')

        setattr(namespace, self.dest, table)

##### Lists #####

class ListFromFile(argparse.Action):
    """Load a list of strings from file"""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ListFromFile, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        with open(values, 'r') as dat:
            flist = dat.readlines()

        setattr(namespace, self.dest, flist)

class ListFromRegex(argparse.Action):
    """Load a list of files from a regular expression"""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ListFromRegex, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from glob import glob
        flist = sorted(glob(os.path.expanduser(values)))

        setattr(namespace, self.dest, flist)

##### Quantities #####

class readQuantity(argparse.Action):
    """Read quantity or a quantity list from the cmd line"""

    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        if nargs < 2 or nargs in ['*', '+', '?']:
            print('WARNING: changed from previous version, nargs are allowed')
        super(readQuantity, self).__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        vals = np.array(values[:-1], dtype=float)
        unit = u.Unit(values[-1])
        if len(vals)==1:
            vals = vals[0]
        vals = vals*unit
        setattr(namespace, self.dest, vals)

##### Advanced processing #####

class PeakPosition(argparse.Action):
    """Load FITS file and get peak postion"""

    def __init__(self, option_strings, dest, nargs='*', **kwargs):
        if nargs not in ['*', '?', '+']:
            raise ValueError('nargs=%s not accepted for PeakToPosition' % nargs)
        super(PeakPosition, self).__init__(option_strings, dest, nargs=nargs,
                **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        positions = []
        values = validate_files(values, check_is_file=True)
        for val in values:
            # Data
            imgi = fits.open(val)[0]
            wcsi = wcs.WCS(imgi.header, naxis=('longitude','latitude'))

            # Maximum
            data = np.squeeze(imgi.data)
            ymax, xmax = np.unravel_index(np.nanargmax(data), data.shape)
            positions += [(xmax, ymax)]

        setattr(namespace, self.dest, positions)

##### Others #####

class NormalizePath(argparse.Action):
    """Normalizes a path or filename"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=False)
        setattr(namespace, self.dest, values)

class CheckFile(argparse.Action):
    """Validates file and check if file or path exist"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values)
        setattr(namespace, self.dest, values)

class startLogger(argparse.Action):
    """Create a logging.logger"""

    def __init__(self, option_strings, dest, nargs=1, **kwargs):
        if nargs!=1:
            raise ValueError("nargs value not allowed")
        default = kwargs.setdefault('default','debug_main.log')
        kwargs['default'] = get_logger('__main__', file_name=default)
        super(startLogger, self).__init__(option_strings, dest, nargs=nargs,
                **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        logger = get_logger('__main__', file_name=values[0])
        setattr(namespace, self.dest, logger)

class readSkyCoords(argparse.Action):
    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        if nargs not in ['*', '+', '?']:
            try:
                if nargs%2!=0:
                    raise ValueError('Only even nargs allowed for coordinates')
            except TypeError:
                raise ValueError('Only even nargs allowed for coordinates')
        kwargs.setdefault('metavar',('RA Dec',)*2)
        super(readSkyCoords, self).__init__(option_strings, dest, nargs=nargs,
                **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values)%2 != 0:
            raise ValueError('Odd number of coordinate values')

        vals = []
        for ra, dec in zip(values[::2], values[1::2]):
            vals += [apycoord.SkyCoord(ra, dec, frame='icrs')]

        setattr(namespace, self.dest, vals)
