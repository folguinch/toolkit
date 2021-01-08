from typing import List, Union
import argparse
import os
import pathlib

from astropy import wcs
from astropy.io import fits
from configparseradv import configparser
import astropy.coordinates as apycoord
import astropy.units as u
import numpy as np
from spectral_cube import SpectralCube

from ..classes.dust import Dust
from ..logger import get_logger

def validate_path(path: pathlib.Path, 
                  check_is_file: bool = False,
                  check_is_dir: bool = False, 
                  mkdir: bool = False):
    """Performs several checks on input path.

    Args:
      filenames: path to check.
      check_is_file: optional; check whether is file and exists.
      check_is_dir: optional; check whther is a directory and exists.
      mkdir: optional; make directories if they do not exist.
    Returns:
      A validated resolved pathlib.Path object
    """
    path = path.expanduser().resolve()
    if check_is_file and not path.is_file():
        raise IOError(f'{path} does not exist')
    elif check_is_dir or mkdir:
        if not path.is_dir() and not mkdir:
            raise IOError(f'{path} directory does not exist')
        else:
            path.mkdir(parents=True, exist_ok=True)

    return path

def validate_paths(filenames: Union[str, List[str]], 
                   check_is_file: bool = False, 
                   check_is_dir: bool = False, 
                   mkdir: bool = False):
    """Performs several checks on input list of file names.

    Args:
      filenames: list of filenames to check.
      check_is_file: optional; check whether is file and exists.
      check_is_dir: optional; check whther is a directory and exists.
      mkdir: optional; make directories if they do not exist.
    Returns:
      Validate path.Path from the input strings.
    """
    try:
        # Case single string file name
        validated = pathlib.Path(filenames)
        validated = validate_path(validated, check_is_file=check_is_file,
                                  check_is_dir=check_is_dir, mkdir=mkdir)
    except TypeError:
        # Case multiple file names
        validated = []
        for filename in filenames:
            aux = pathlib.Path(filename)
            aux = validate_path(aux, check_is_file=check_is_file,
                                check_is_dir=check_is_dir, mkdir=mkdir)
            validated.append(aux)

    return validated

# Loader actions
class LoadConfig(argparse.Action):
    """Action class for loading a configuration file in argparse."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=True)
        config = configparser.ConfigParserAdv()
        aux = config.read(values)
        setattr(namespace, self.dest, config)

class LoadArray(argparse.Action):
    """Action class for loading a np.array from command line."""

    def __call__(self, parser, namespace, values, option_string=None):
        array = np.array(values, dtype=float)
        setattr(namespace, self.dest, array)

class ArrayFromRange(argparse.Action):
    """Action class for creating a np.array with linspace from command line"""

    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        if nargs not in range(2, 5):
            raise ValueError('only 2, 3 or 4 nargs allowed')
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        start, stop = float(values[0]), float(values[1])
        if len(values) == 4:
            n = int(values[2])
            if values[-1] == 'linear':
                fn = np.linspace
            elif values[-1] == 'log':
                fn = np.logspace
                start = np.log10(start)
                stop = np.log10(stop)
            else:
                raise NotImplementedError(f'{values[-1]} not implemented')
        elif len(values) == 3:
            fn = np.linspace
            n = int(values[2])
        else:
            fn = np.linspace
        value = fn(start, stop, n)
        setattr(namespace, self.dest, value)

class LoadStructArray(argparse.Action):
    """Load an structured np.array from file."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from ..array_utils import load_struct_array
        array = load_struct_array(validate_files(values, check_is_file=True))
        setattr(namespace, self.dest, array)

class LoadMixedStructArray(argparse.Action):
    """Load a mixed structured np.array from file."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from ..array_utils import load_mixed_struct_array
        array = load_mixed_struct_array(validate_files(values,
                                                       check_is_file=True))
        setattr(namespace, self.dest, array)

class LoadTXTArray(argparse.Action):
    """Load an np.array from file."""

    def __call__(self, parser, namespace, values, option_string=None):
        array = np.loadtxt(validate_files(values, check_is_file=True), 
                           dtype=float)
        setattr(namespace, self.dest, array)

class LoadFITS(argparse.Action):
    """Action for loading a FITS file with astropy"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=True)
        try:
            vals = fits.open(values.resolve())[0]
        except AttributeError:
            vals = []
            for val in values:
                vals += [fits.open(val)[0]]
        setattr(namespace, self.dest, vals)

class LoadCube(argparse.Action):
    """Action for loading an SpectralCube"""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=True)
        try:
            vals = SpectralCube.read(values.resolve())
        except AttributeError:
            vals = []
            for val in values:
                vals += [SpectralCube.read(val)]
        setattr(namespace, self.dest, vals)

class LoadDust(argparse.Action):
    """Action for loading dust files"""

    def __call__(self, parser, namespace, values, option_string=None):
        dust = Dust(values)
        setattr(namespace, self.dest, dust)

class LoadTable(argparse.Action):
    """Action for loading astropy Tables"""

    def __call__(self, parser, namespace, values, option_string=None):
        from ..tables import Table

        try:
            tabname = ''+values
            #table_id = os.path.splitext(os.path.basename(tabname))[0]
            table = Table(tabname)
        except TypeError:
            if len(values) == 2:
                table = Table(values[0], table_id=values[1])
            elif len(values) == 1:
                tabname = values[0]
                #table_id = os.path.splitext(os.path.basename(tabname))[0]
                table = Table(tabname)
            else:
                raise ValueError('Number of values not allowed.')

        setattr(namespace, self.dest, table)

# Lists actions
class ListFromFile(argparse.Action):
    """Load a list of strings from file."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ListFromFile, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        with open(values, 'r') as dat:
            flist = dat.readlines()

        setattr(namespace, self.dest, flist)

class ListFromRegex(argparse.Action):
    """Load a list of files from a regular expression."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ListFromRegex, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from glob import glob
        flist = sorted(glob(os.path.expanduser(values)))

        setattr(namespace, self.dest, flist)

# Quantity actions
class readQuantity(argparse.Action):
    """Read quantity or a quantity list from the cmd line."""

    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        try:
            if nargs < 2:
                raise ValueError('nargs cannot be < 2')
        except TypeError:
            pass

        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) < 2:
            raise ValueError(f'Cannot read quantity from values: {values}')
        vals = np.array(values[:-1], dtype=float)
        unit = u.Unit(values[-1])
        if len(vals) == 1:
            vals = vals[0]
        vals = vals * unit
        setattr(namespace, self.dest, vals)

class readUnit(argparse.Action):
    """Read quantity or a quantity list from the cmd line."""

    def __call__(self, parser, namespace, values, option_string=None):
        vals = []
        try:
            for val in values:
                vals.append(u.Unit(val))
            if len(vals) == 1:
                vals = vals[0]
        except ValueError:
            vals = u.Unit(values)
        setattr(namespace, self.dest, vals)

# Advanced processing actions
class PeakPosition(argparse.Action):
    """Load FITS file and get peak postion."""

    def __init__(self, option_strings, dest, nargs='*', **kwargs):
        if nargs not in ['*', '?', '+']:
            raise ValueError(f'nargs={nargs} not accepted for PeakPosition')
        kwargs.setdefault('metavar', ('FITSFILE',)*2)
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        positions = []
        values = validate_files(values, check_is_file=True)
        for val in values:
            # Data
            imgi = fits.open(val)[0]

            # Maximum
            data = np.squeeze(imgi.data)
            ymax, xmax = np.unravel_index(np.nanargmax(data), data.shape)
            positions += [(xmax, ymax)]

        setattr(namespace, self.dest, positions)

class readSkyCoords(argparse.Action):
    """Read one or more sky coordinates."""

    def __init__(self, option_strings, dest, nargs=2, **kwargs):
        try:
            if nargs < 2:
                raise ValueError('Only nargs values >= 2 accepted')
            if nargs%2 == 0:
                kwargs.setdefault('metavar', ('RA Dec',)*nargs)
            else:
                kwargs.setdefault('metavar', 
                                  ('RA Dec ',)*(nargs-1) + ('FRAME',))
        except TypeError:
            kwargs.setdefault('metavar', ('RA Dec',)*2 + ('[FRAME]',))
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) < 2:
            raise ValueError('Could not read sky coordinate')
        elif len(values)%2 == 0:
            frame = 'icrs'
        else:
            frame = values[-1]
            values = values[:-1]
        vals = []
        for ra, dec in zip(values[::2], values[1::2]):
            vals.append(apycoord.SkyCoord(ra, dec, frame=frame))

        setattr(namespace, self.dest, vals)

# Path actions
class NormalizePath(argparse.Action):
    """Normalizes a path or filename."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values)
        setattr(namespace, self.dest, values)

class MakePath(argparse.Action):
    """Check and create directory if needed."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, mkdir=True)
        setattr(namespace, self.dest, values)

class CheckFile(argparse.Action):
    """Validates files and check if they exist."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = validate_files(values, check_is_file=True)
        setattr(namespace, self.dest, values)

# Loger actions
class startLogger(argparse.Action):
    """Create a logger."""

    def __init__(self, option_strings, dest, nargs=1, **kwargs):
        if nargs not in [0, 1]:
            raise ValueError("nargs value not allowed")
        default = kwargs.setdefault('default', 'debug_main.log')
        if nargs == 0:
            kwargs['default'] = get_logger('__main__', filename=default)
        else:
            kargs.setdefault('metavar', 'LOGFILE')
        # Testing
        print(option_strings)
        super().__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 1:
            logger = get_logger('__main__', filename=values[0])
        setattr(namespace, self.dest, logger)

