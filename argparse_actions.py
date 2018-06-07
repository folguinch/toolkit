import os, argparse
from configparser import ExtendedInterpolation

import numpy as np
from astropy.io import fits

from .myconfigparser import myConfigParser
from .classes.dust import Dust

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

class LoadStructArray(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(LoadStructArray, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from .array_utils import load_struct_array
        array = load_struct_array(validate_files(values))
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
        flist = glob(values)

        setattr(namespace, self.dest, flist)

##### Others #####

class NormalizePath(argparse.Action):
    """Normalizes a path or filename"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=False)
        setattr(namespace, self.dest, values)

class CheckFile(argparse.Action):
    """Normalizes a path or filename"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values)
        setattr(namespace, self.dest, values)
