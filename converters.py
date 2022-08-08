"""Convert between data types."""
from typing import Dict, Union, Optional
import argparse
import configparser

from astropy import wcs
from astropy.io import fits
import astropy.units as u
import numpy as np

def argparser_to_configparser(
    args: argparse.Namespace,
    config: configparser.ConfigParser,
    section: str,
) -> configparser.ConfigParser:
    """Update a `configparser` object from the values in `args`.

    If option is in `args` but not in `config`, then it is skipped.

    Args:
      args: argument parser object.
      config: configuration parser.
      section: section of `config` to update

    Returns:
      An updated configuration parser.
    """
    # Initial values
    new_values = {}
    args_dict = vars(args)

    # Filter values
    for key, val in args_dict.items():
        if key in config.options(section):
            new_values[key] = val

    # Update config parser
    config.read_dict({section: new_values})

    return config

def array_to_hdu(array: Union[np.array, u.Quantity],
                 reference: Union[fits.PrimaryHDU, Dict],
                 unit: Optional[u.Unit] = None) -> fits.PrimaryHDU:
    """Convert a `np.array` to a `PrimaryHDU`.

    If `array` is a `Quantity` array and `unit` is given, then `array` is
    converted to unit. Otherwise, if `unit` is given then it assumed to be the
    unit of `array`.

    Args:
      array: data to convert.
      reference: reference image to obtain WCS information.
      unit: optional; unit or output unit of the data.
    """
    # Get header
    if hasattr(reference, 'wcs'):
        header = reference.wcs.to_header()
    else:
        header = wcs.WCS(reference, naxis=['longitude', 'latitude']).to_header()

    # Convert data
    if hasattr(array, 'unit'):
        if unit is None:
            data = np.squeeze(array.value)
            unit = array.unit
        else:
            data = np.squeeze(array.to(unit).value)
        header['BUNIT'] = f'{unit:FITS}'
    else:
        if unit is not None:
            header['BUNIT'] = f'{unit:FITS}'
        data = np.squeeze(array)

    return fits.PrimaryHDU(data=data, header=header)

def quantity_from_hdu(hdu: fits.PrimaryHDU,
                      unit: Optional[u.Unit] = None) -> u.Quantity:
    """Convert a `PrimaryHDU` to a `Quantity`."""
    # Data
    data = np.squeeze(hdu.data)

    # Unit
    bunit = u.Unit(hdu.header['BUNIT'])
    data = data * bunit
    if unit is not None:
        data = data.to(unit)

    return data

