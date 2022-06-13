"""Functions for working with numpy arrays."""
from typing import List, Optional, Tuple, Union
import pathlib

import numpy as np
import astropy.units as u

def check_composed_units(line: List[str],
                         left_delimiter: str = '[',
                         right_delimiter: str = ']') -> List[str]:
    """Check composed units from strings in file headers.

    Args:
      line: text line splitted.
      left_delimiter: optional; left char delimiter for composed unit.
      right_delimiter: optional; right char delimiter for composed unit.

    Return:
      Validated units.
    """
    # Check delimiters
    if len(left_delimiter) != 1 or len(right_delimiter) != 1:
        raise ValueError('Delimiter length must be 1')

    # Process units
    units = []
    aux = ''
    for unit in line:
        if left_delimiter in unit or aux != '':
            aux = f'{aux} {unit}'
            if aux[-1] == right_delimiter:
                aux = aux[1:-1]
                units.append(aux)
                aux = ''
            else:
                aux = aux.strip()
        else:
            units.append(unit)

    return units

def _load_struct_array(
    file_name: pathlib.Path,
    usecols: Optional[List[int]] = None,
    nounit: str = '-',
    empty_unit: Optional[u.Unit] = None,
    dtype: Optional[np.dtype] = None,
) -> Tuple[np.array, dict]:
    """Load a structured array table.

    Each file has a header. The first row has the name of each parameter and
    the second the units. The name of each parameter must be different. String
    values or dimensionless are marked with unit `-` in the second line.

    If dtype is given, np.loadtxt is used else the dtype are determined by
    np.genfromtxt.

    Args:
      file_name: file to be loaded.
      usecols: optional; columns to load.
      nounit: optional; char for columns without units.
      empty_unit: optional; unit for nounit column.
      dtype: optional; data type.

    Returns:
      A tuple composed of:

        - data: structured array.
        - units: dictionary with the units.

    Raises:
      ValueError if header names have spaces or different number of column
        names and number of columns.
    """
    with file_name.open('r') as inp:
        # Read first two lines
        line1 = inp.readline().strip(' #').split()
        line2 = inp.readline().strip(' #').split()
        line2 = check_composed_units(line2)

        # Check header
        if len(line1) != len(set(line1)):
            raise ValueError('Header length different than number of columns')

        # Read dtype and data units
        names = []
        units = {}
        for i, (name, unit) in enumerate(zip(line1, line2)):
            # Skip data column
            if usecols is not None and i not in usecols:
                continue

            # Names
            if dtype is None:
                names.append(name)
            else:
                names.append((name, dtype))

            # Units
            if unit == nounit:
                units[name] = empty_unit
            else:
                units[name] = u.Unit(unit)

        # Read remaininf data data
        if dtype is None:
            data = np.genfromtxt(inp, usecols=usecols, dtype=None, names=names)
        else:
            data = np.loadtxt(inp, usecols=usecols, dtype=names)

    return data, units

def filename_to_list(
    file_name: Union[pathlib.Path, List[pathlib.Path]],
) -> Union[list, bool]:
    """Determine if there are multiple file names."""
    if not hasattr(file_name, 'expanduser'):
        file_names = file_name
        was_list = True
    else:
        file_names = [file_name]
        was_list = False

    return file_names, was_list

def load_struct_array(
    file_name: Union[pathlib.Path, List[pathlib.Path]],
    usecols: Optional[List[int]] = None,
    nounit: str = '-',
    empty_unit: Optional[u.Unit] = u.Unit(1),
) -> Union[Tuple[np.array, dict], List[Tuple[np.array, dict]]]:
    """Load a structured array table.

    Args:
      file_name: file or files to be loaded.
      usecols: optional; columns to load.
      nounit: optional; char for columns without units.
      empty_unit: optional; unit for nounit column.
    Returns:
      data: structured array.
      units: dictionary with the units.
    """
    data = []
    files, was_list = filename_to_list(file_name)
    for f in files:
        aux = _load_struct_array(f, usecols=usecols, nounit=nounit,
                                 empty_unit=empty_unit, dtype=float)
        data.append(aux)

    if not was_list:
        return data[0]
    else:
        return data

def load_mixed_struct_array(
    file_name: Union[pathlib.Path, List[pathlib.Path]],
    usecols: Optional[List[int]] = None,
    nounit: str = '-',
    empty_unit: Optional[u.Unit] = None,
) -> Union[Tuple[np.array, dict], List[Tuple[np.array, dict]]]:
    """Load a structured array table of mixed types.

    Args:
      file_name: file or files to be loaded.
      usecols: optional; columns to load.
      nounit: optional; char for columns without units.
      empty_unit: optional; unit for nounit column.
    Returns:
      data: structured array.
      units: dictionary with the units.
    """
    data = []
    files, was_list = filename_to_list(file_name)
    for f in files:
        aux = _load_struct_array(f, usecols=usecols, nounit=nounit,
                                 empty_unit=empty_unit)
        data.append(aux)

    if not was_list:
        return data[0]
    else:
        return data

def save_struct_array(file_name: pathlib.Path, data: np.array, units: dict,
                      fmt: str = '%10.4e\t') -> None:
    """Save a structured array table.

    Save the data in a way it can be loaded by the load_struct_array function.

    Args:
      file_name: name of the file.
      data: array to save.
      units: physical units of the data.
      fmt: optional; string format for the data.
    """
    # Lines
    line1 = []
    line2 = []
    for name in data.dtype.names:
        line1.append(f'{name}')
        if units[name] is None or units[name] == u.Unit(''):
            unit = '-'
        else:
            unit = f'{units[name]}'
        if ' ' in unit:
            unit = f'[{unit}]'
        line2.append(f'{unit}')
    line1 = '#' + '\t'.join(line1).strip() + '\n'
    line2 = '#' + '\t'.join(line2).strip() + '\n'
    lines = line1 + line2

    # Data
    try:
        lines += '\n'.join((fmt * len(d)).strip() % tuple(d) for d in data)
    except TypeError:
        lines += '\n'.join(fmt.strip() % tuple(d) for d in data)
    file_name.write_text(lines)

