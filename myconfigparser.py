import os
from configparser import ConfigParser

import numpy as np
import astropy.units as u

class myConfigParser(ConfigParser):
    """Extend the configparser.ConfigParser behaviour.
    """

    def getquantity(self, *args, **kwargs):
        """Return an ``astropy.Quantity`` from the parser data.

        Input parameters are the same as for ``parser.get(*args,**kwargs)``
        """
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            val= val.split()
        if len(val)==2:
            return float(val[0]) * u.Unit(val[1])
        else:
            return np.array(val[:-1], dtype=float) * u.Unit(val[-1])

    def getlist(self, *args, **kwargs):
        """Return a list of strings"""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            if ' ' not in val:
                return [val]
            else:
                return val.split()

    def getfloatlist(self, *args, **kwargs):
        """Return a list of float values."""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            if ' ' not in val:
                val = [val]
            else:
                val= val.split()
        return map(float, val)

    def getpath(self, *args, **kwargs):
        """Return a real path"""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            return os.path.expanduser(os.path.expandvars(val))
