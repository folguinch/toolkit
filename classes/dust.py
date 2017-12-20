import astropy.units as u

from .data_1d import Data1D
from ..logger import get_logger

class Dust(Data1D):
    """Create a Dust objects.

    Attributes:
        address: file name.
        data: the data.
        units: units of the elements of the data.
        logger: logging manager.
    """
    logger = get_logger(__name__)

    def interpolate(self, newx, **kwargs):
        """Interpolate and evaluate in *newx*.

        Parameters:
            newx (astropy.quantity): the new wavelength values.
            kwargs: parameters for scipy.interp1d
        """
        cnewx = newx.to(self.units['wlg'], equivalencies=u.spectral())
        return super(Dust, self).interpolate('wlg', 'kappa', cnewx, **kwargs)

