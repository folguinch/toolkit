import os, argparse
from configparser import ExtendedInterpolation

import numpy as np

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
            validated += [os.path.expanduser(filenames)]
            if check_is_file and not os.path.isfile(validated[-1]):
                raise IOError('%s does not exist' % (validated[-1]))

    return validated

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
