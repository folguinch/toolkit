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

        # If the value was changed then it may be already a quantity
        if hasattr(val, 'unit'):
            return val
        
        # Fallback values
        if val is None or len(val)==0:
            return val
        else:
            val = val.split()

        # Convert string to quantity
        if len(val)==1:
            # Dimesionless
            return float(val[0])
        elif len(val)==2:
            # Single quantity
            return float(val[0]) * u.Unit(val[1])
        else:
            # Array of values
            return np.array(val[:-1], dtype=float) * u.Unit(val[-1])

    def getlist(self, *args, **kwargs):
        """Return a list of strings"""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            try:
                return val.split()
            except AttributeError:
                return list(val)
            #if ' ' not in val:
            #    return [val]
            #else:
            #    return val.split()

    def getfloatlist(self, *args, **kwargs):
        """Return a list of float values."""
        #val = self.get(*args, **kwargs)
        #if val is None:
        #    return val
        #else:
        #    if ' ' not in val:
        #        val = [val]
        #    else:
        #        val= val.split()
        val = self.getlist(*args, **kwargs)
        if val is None:
            return val
        else:
            return map(float, val)

    def getintlist(self, *args, **kwargs):
        """Return a list of int values."""
        #val = self.get(*args, **kwargs)
        #if val is None:
        #    return val
        #else:
        #    if ' ' not in val:
        #        val = [val]
        #    else:
        #        val= val.split()
        val = self.getlist(*args, **kwargs)
        if val is None:
            return val
        else:
            return map(int, val)

    def getpath(self, *args, **kwargs):
        """Return a real path"""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            return os.path.expanduser(os.path.expandvars(val))
