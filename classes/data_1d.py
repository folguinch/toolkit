from scipy.interpolate import interp1d

from ..array_utils import *
from ..logger import get_logger
from .data import Data

class Data1D(Data):
    """Creates a 1-D data structure.

    This object is useful to managa data like radial profiles, spectra or
    slices.

    Attributes:
        address: file name.
        data: the data.
        units: units of the elements of the data.
        logger: logging manager.
    """

    def __init__(self, file_name=None, data=None, units=None, wlg=None):
        """Creates a new profile object.

        Parameters:
            file_name (str): file name of the profile.
            wlg (float, default=None): wavelength.
        """
        self.units = units
        super(Data1D, self).__init__(file_name, data)
        self.logger = get_logger(__name__)

    def __getitem__(self, key):
        return self.data[key] * self.units[key]

    def __setitem__(self, key, val):
        assert hasattr(val, 'value')
        assert hasattr(val, 'unit')
        assert len(val.value)==len(self[key].value)
        self.data[key] = val.value
        self.units[key] = 1. * val.unit

    def load(self):
        """Load data from file.

        Each file has a standard structure with the data in different columns.
        """
        self.data, self.units = load_struct_array(self.address)

    def save(self, file_name=None):
        """Save the data to file.

        Parameters:
            file_name (default=None): new file name.
        """
        save_struct_array(file_name, self.data, self.units)

    def as_function(self, keyx, keyy, **kwargs):
        """Interpolate the data and return a function.

        Parameters:
            keyx: key of the x value.
            keyy: key of the y value.
            kwargs: arguments for scipy.interpolate.interp1d
        """
        kind = kwargs.pop('kind', 'linear')
        bound_error = kwargs.pop('bounds_error', False)
        fill_value = kwargs.pop('fill_value', 0.)
        return interp1d(self.data[keyx], self.data[keyy], kind=kind,
                bounds_error=bound_error, fill_value=fill_value, **kwargs)

    def interpolate(self, keyx, keyy, newx, **kwargs):
        """Interpolate the data and evaluate in new points.

        Parameters:
            keyx: key of the x data.
            keyy: key of the y data.
            newx: values where the interpolated function is evaluated.
            kwargs: arguments for scipy.interpolate.interp1d
        """
        fn = self.as_function(keyx, keyy, **kwargs)
        return fn(newx.to(self.units[keyx]).value) * self.units[keyy]

    def convert(self, key, new_unit):
        """Convert the units of data in *key* to *new_unit*.

        Parameters:
            key (str): data to convert.
            new_unit (astropy.unit): new physical unit.
        """
        self.data[key] = (self.data[key]*self.units[key]).to(new_unit).value
        self.units[key] = 1.*new_unit
