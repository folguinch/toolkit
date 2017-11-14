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
        val = self.get(*args, **kwargs).split()
        if len(val)==2:
            return float(val[0]) * u.Unit(val[1])
        else:
            return np.array(val[:-1], dtype=float) * u.Unit(val[-1])
