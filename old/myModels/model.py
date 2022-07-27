from collections import OrderedDict
from abc import ABCMeta, abstractmethod

import numpy as np
from astropy.modeling.fitting import LevMarLSQFitter, _convert_input, _model_to_fit_params

class Model(object):
    """Model ABC.

    Attributes:
        params (dict): model parameters.
        units (dict): standard units for the model.
        constraints (tuple): constraints on model parameters.
    """
    __metaclass__ = ABCMeta

    constraints = (-np.inf, np.inf)

    def __init__(self, params):
        """Initial parameter of the model.

        The units of the initial parameters are the same for all the values
        used in the fitting routines.
        
        Parameters:
            params (OrderedDict): parameters of the model with units.
        """
        self.validate(params)
        self.params = OrderedDict()
        self.units = {}
        self._fill_values(params)

    def __getitem__(self, key):
        return self.params[key] * self.units.get(key, 1.)

    def __setitem__(self, key, val):
        self._fill_values({key:val}, check_existance=True)

    def __call__(self, x, *args, **kwargs):
        self._update_values(*args)
        return self.model_function(x, **kwargs)

    @abstractmethod
    def validate(self, params):
        pass

    @abstractmethod
    def model_function(self, *args, **kwargs):
        pass

    @property
    def values(self):
        return self.params.values()

    @property
    def bounds(self):
        try:
            bounds = np.array(self.constraints.values(),
                    dtype=[('lower',float),('upper',float)])
            return bounds['lower'], bounds['upper']
        except AttributeError:
            if len(self.constraints)==2:
                return self.constraints

    def _fill_values(self, params, check_existance=False):
        """Fill the attributes."""
        for key,val in params.items():
            if check_existance and key not in self.params.keys():
                raise KeyError('This model does not have parameter %s' % key)

            if hasattr(val, 'unit'):
                self.params[key] = val.value
                self.units[key] = val.unit
            else:
                self.params[key] = val

    def _update_values(self, *args):
        """Update all the parameter values."""
        if len(args) == 0:
            pass
        else:
            assert len(args)==len(self.params)
            assert not all([hasattr(p, 'unit') for p in args])
            self.params = OrderedDict(zip(self.params.keys(),args))


