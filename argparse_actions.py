import os, argparse
from configparser import ExtendedInterpolation

import numpy as np

from .myconfigparser import myConfigParser

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
