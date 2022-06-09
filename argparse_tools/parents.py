"""Argparse parent parsers commonly used."""
from typing import Callable, Optional, Union
import argparse
import pathlib
import warnings

try:
    import astro_source.source as source
except ImportError:
    print('astro_source not available')
    pass

from .functions import positions_to_pixels
import toolkit.argparse_tools.actions as actions

# Typing
PosFunction = Callable[[argparse.Namespace], None]
Path = Union[pathlib.Path, str]

def astro_source(parser: argparse.ArgumentParser) -> None:
    """Read an `astro_source.Source`."""
    try:
        parser.add_argument('--source',
                            action=source.LoadSources,
                            help='Source(s) configuration file(s)')
    except NameError:
        pass

def source_position(
    required: bool = False,
    function: PosFunction = positions_to_pixels,
) -> argparse.ArgumentParser:
    """Parent for reading a source position.

    The function creates 2 defaults:

      - `pos`: a list of position pairs. The values stored ar in xy pixels.
      - `position_fn`: which generates and fills the `pos` from the command
        line input.

    Args:
      required: optional; is the argument required?
      function: optional; a function mapping the postion values to xy pixels.
    """
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
    warnings.warn(('The values pos and position_fn will be removed from'
                   'the parser in future versions of this parent'),
                  FutureWarning)
    parser.set_defaults(position_fn=function, pos=[])

    return parser

def logger(filename: Optional[Path] = None) -> argparse.ArgumentParser:
    """Parent parser to initiate a logging system.

    Args:
      filename: optional; default filename for logging.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-v', '--vv', '--vvv', default=filename,
                        action=actions.StartLogger,
                        help='Logging setup')

    return parser

def verify_files(*args, **kwargs) -> argparse.ArgumentParser:
    """Create a parser with the input strings that verify for file
    existance.

    Args:
      args: argument name.
      kwargs: additional arguments for `add_argument` for each `arg`.
    """
    parser = argparse.ArgumentParser(add_help=False)
    for opt in args:
        if opt.startswith('-'):
            key = opt.strip('-')
        else:
            key = opt
        parser.add_argument(f'{opt}', action=actions.CheckFile,
                            **kwargs[key])

    return parser

def paths(*args, **kwargs) -> argparse.ArgumentParser:
    """Create a parser with the input strings that create paths.

    Args:
      args: argument name.
      kwargs: additional arguments for `add_argument` for each `arg`.
    """
    parser = argparse.ArgumentParser(add_help=False)
    for opt in args:
        parser.add_argument(f'--{opt}', action=actions.MakePath,
                            **kwargs[opt])

    return parser
