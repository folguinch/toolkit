import os
from configparser import ConfigParser

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

def converters(value, dtype):
    if dtype.lower()=='quantity':
        aux = value.split()
        if len(aux) == 2:
            return float(aux[0])*u.Unit(aux[1])
        else:
            return np.array(aux[:-1], dtype=float) * u.Unit(aux[-1])
    elif dtype.lower()=='skycoord':
        try:
            ra, dec, frame = value.split()
        except ValueError:
            ra, dec = value.split()
            frame = 'icrs'
        return SkyCoord(ra, dec, frame=frame)
    else:
        raise NotImplementedError('converter to %s not available' % dtype)
 
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
        else:
            # Or it may be a dimensionless float or float array
            try:
                # Set dimenssionless unit
                return (val + 0) * u.Unit(1)
            except TypeError:
                pass
        
        # Fallback values
        if val is None or len(val)==0:
            return val
        else:
            val = val.split()

        # Convert string to quantity
        if len(val)==1:
            # Dimesionless
            return float(val[0]) * u.Unit(1)
        elif len(val)==2:
            # Single quantity
            return float(val[0]) * u.Unit(val[1])
        else:
            # Array of values
            try:
                # Check if dimensionless
                aux = float(val[-1])
                return np.array(val, dtype=float) * u.Unit(1)
            except ValueError:
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

    def getskycoord(self, *args, **kwargs):
        """Return an astropy SkyCoord"""
        val = self.get(*args, **kwargs)
        if val is None:
            return val
        else:
            ra, dec, frame = val.split()
            return SkyCoord(ra, dec, frame=frame)

    def getvalue(self, *args, **kwargs):
        """Get values"""
        # Options
        opts = {'fallback':None, 'n':None, 'sep':' ', 'dtype':None,
                'allow_global':True}
        opts.update(kwargs)
        
        # Abbreviations
        n = opts['n']
        key = args[1]
        getkw = {'fallback':opts['fallback']}

        # Get values
        if n is not None and (key + str(n)) in self.options(args[0]):
            newkey = key + str(n)
            value = self.get(args[0], newkey, **getkw)
        elif key in self.options(args[0]):
            value = self.get(*args, **getkw)
            if n is not None:
                value = value.split(opts['sep'])
                if len(value)==1:
                    if n!=0 and not opts['allow_global']:
                        value = opts['fallback']
                    else:
                        value = value[0]
                else:
                    try:
                        value = value[n]
                    except IndexError:
                        print('WARNING: %s not in values list, using fallback'\
                        % key)
                        value = opts['fallback']
        else:
            value = opts['fallback']

        try:
            value = value.strip()
        except:
            pass

        if opts['dtype'] is None or value is None:
            return value
        else:
            try:
                return opts['dtype'](value)
            except TypeError:
                return converters(value, opts['dtype'])

    def getvalueiter(self, *args, **kwargs):
        opts = {'sep':kwargs.get('sep', ' '), 'n':0, 'allow_global':False}
        while self.getvalue(*args, **opts):
            kwargs['n'] = opts['n']
            value = self.getvalue(*args, **kwargs)
            if value is not None:
                yield value
            else:
                break
            opts['n'] = opts['n'] + 1


