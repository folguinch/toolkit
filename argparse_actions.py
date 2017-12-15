import os, argparse
from configparser import ExtendedInterpolation

import numpy as np
from astroSource.image import Image

from .myconfigparser import myConfigParser
from .dust import Dust

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
    """Action class for loading fits files with mySoRadfit Image"""

    def __call__(self, parser, namespace, values, option_string=None):
        images = []
        try:
            assert os.path.isfile(values)
            images += [Image(values)]
        except TypeError:
            for val in values:
                assert os.path.isfile(val)
                images += [Image(val)]
        setattr(namespace, self.dest, images)

class LoadDust(argparse.Action):
    """Action for loading dust files"""

    def __call__(self, parser, namespace, values, option_string=None):
        dust = Dust(values)
        setattr(namespace, self.dest, dust)
