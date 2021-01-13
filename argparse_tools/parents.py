from typing import Callable, Optional, TypeVar
import argparse

import actions
import functions as fns

"""Argparse parent parsers commonly used"""

# Typing
PosFunction = Callable[[argparse.Namespace], None]
Path = TypeVar('Path')

def astro_source(parser):
    try:
        import astroSource.source as source
        parser.add_argument('--source', action=source.LoadSourcesfromConfig,
                help='Source(s) configuration file(s)')
    except ImportError:
        print('astroSource not available')
        pass

def source_position(required: bool = False, 
                    function: PosFunction = fns.positions_to_pos):

    parser = argparse.ArgumentParser(add_help=False)
    group1 = parser.add_mutually_exclusive_group(required=required)
    group1.add_argument('--coordinate', nargs='*', 
            action=actions.ReadSkyCoords,
            help='Sky coordinates with units or : separator')
    group1.add_argument('--position', metavar=('X Y',)*2, nargs='*', type=int,
            help='Positions')
    group1.add_argument('--reference', metavar='IMG',
            action=actions.PeakPosition,
            help='Reference image to get position from peak')
    astro_source(group1)
    parser.set_defaults(position_fn=function, pos=[])

    return parser

def logger(filename: Optional[str, Path] = None) -> argparse.ArgumentParser:
    """Parent parser to initiate a logging system.

    Args:
      filename: optional; default filename for logging.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-v', '--vv', '--vvv', 
                        default=filename,
                        action=actions.StartLogger,
                        help='Logging setup')

    return parser
