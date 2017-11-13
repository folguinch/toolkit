import argparse
from configparser import ConfigParser, ExtendedInterpolation

class LoadConfig(argparse.Action):
    """Action class for loading a configuration file in argparse"""

    def __call__(self, parser, namespace, values, option_string=None):
        config = ConfigParser(interpolation=ExtendedInterpolation())
        config.read(values)
        print 'Read config %s' % values
        setattr(namespace, self.dest, config)

