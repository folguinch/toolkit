"""Argparse parent parsers commonly used."""
from typing import Callable, Optional, Union
import argparse
import pathlib
import warnings

try:
    from astro_source import source
except ImportError:
    print('astro_source not available')

from .functions import positions_to_pixels
from .actions import (ReadSkyCoords, PeakPosition, StartLogger, CheckFile,
                      MakePath, ReadQuantity)

# Typing
PosFunction = Callable[[argparse.Namespace], None]
Path = Union[pathlib.Path, str]

def _astro_source(parser: argparse.ArgumentParser) -> None:
    """Read an `astro_source.Source`."""
    try:
        parser.add_argument('--source',
                            action=source.LoadSources,
                            help='Source(s) configuration file(s)')
    except NameError:
        pass

def astro_source() -> argparse.ArgumentParser:
    """Read an `source.Source` and return it in parser."""
    parser = argparse.ArgumentParser(add_help=False)
    _astro_source(parser)
    return parser

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
      required: Optional. Is the argument required?
      function: Optional. A function mapping the postion values to xy pixels.
    """
    parser = argparse.ArgumentParser(add_help=False)
    group1 = parser.add_mutually_exclusive_group(required=required)
    group1.add_argument('--coordinate', nargs='*',
                        action=ReadSkyCoords,
                        help='Sky coordinates with units or : separator')
    group1.add_argument('--position', metavar=('X Y',)*2, nargs='*', type=int,
                        help='Positions')
    group1.add_argument('--reference', metavar='IMG',
                        action=PeakPosition,
                        help='Reference image to get position from peak')
    _astro_source(group1)
    #warnings.warn(('The values pos and position_fn will be removed from'
    #               'the parser in future versions of this parent'),
    #              FutureWarning)
    parser.set_defaults(position_fn=function)

    return parser

def source_properties(properties: Sequence[str],
                      add_source: bool = False,
                      required: bool = False) -> argparse.ArgumentParser:
    """Parent for source properties as quantities.

    Args:
      properties: Source properties requested.
      add_source: Optional. Include `source.Source`?
      required: Optional. Are all properties required?
    """
    parser = argparse.ArgumentParser(add_help=False)
    if add_source:
        _astro_source(parser)
    for prop in properties:
        if required:
            opt = f'{prop}'
        else:
            opt = f'--{prop}'
        parser.add_argument(opt, metavar=('VAL', 'UNIT'), default=None,
                            action=ReadQuantity,
                            help=f'Source {prop}')

    return parser

def logger(filename: Optional[Path] = None) -> argparse.ArgumentParser:
    """Parent parser to initiate a logging system.

    Args:
      filename: optional; default filename for logging.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-v', '--vv', '--vvv', default=filename,
                        action=StartLogger,
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
        parser.add_argument(f'{opt}', action=CheckFile,
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
        parser.add_argument(f'--{opt}', action=MakePath,
                            **kwargs[opt])

    return parser
