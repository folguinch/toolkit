import argparse

import actions
import functions as fns
from ..logger import get_logger

"""Argparse parent parsers commonly used"""

def astro_source(parser):
    try:
        import astroSource.source as source
        parser.add_argument('--source', action=source.LoadSourcesfromConfig,
                help='Source(s) configuration file(s)')
    except ImportError:
        print('astroSource not available')
        pass

def source_position(required=False, function=fns.positions_to_pos):

    parser = argparse.ArgumentParser(add_help=False)
    group1 = parser.add_mutually_exclusive_group(required=required)
    group1.add_argument('--coordinate', nargs='*', 
            action=actions.readSkyCoords,
            help='Sky coordinates with units or : separator')
    group1.add_argument('--position', metavar=('X Y',)*2, nargs='*', type=int,
            help='Positions')
    group1.add_argument('--reference', metavar='IMG',
            action=actions.PeakPosition,
            help='Reference image to get position from peak')
    astro_source(group1)
    parser.set_defaults(position_fn=function, pos=[])

    return parser

def logger():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--loglevel', default=['info'], nargs=1,
            choices=['info','debug','warn','error'],
            help='Logging stdout level')
    parser.set_defaults(log=get_logger)

    return parser
