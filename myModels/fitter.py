from scipy.optimize import least_squares
import numpy as np

class Fitter(object):
    """Manage fitting of models.

    Attributes:
        results: fitting results.
    """

    def __init__(self):
        self.results = None

    def fit_lsq(self, model, x, y, z=None, err=None, **kwargs):
        """Fit the model to the given data.

        Parameters:
            model: model object.
            x: coordinate data.
            y: coordinate or observed data.
            z: observed data.
            err: observed data errors.
        """
        self.results = least_squares(self.fit_func, model.values, 
                bounds=model.bounds, args=(model,),
                kwargs={'x':x, 'y':y, 'z':z, 'err':err},
                jac=self.deriv_func, **kwargs)

    @staticmethod
    def fit_func(params, *args, **kwargs):
        """Function to minimize."""
        model = args[0]
        x = kwargs['x']
        y = kwargs['y']
        z = kwargs.get('z')
        if z is not None:
            err = kwargs.get('err', 1.*z.unit).to(z.unit)
            mod = model(x, y, *params).to(z.unit)
            dat = z
        else:
            err = kwargs.get('err', 1.*y.unit).to(y.unit)
            mod = model(x, *params).to(y.unit)
            dat = y

        return np.ravel((dat.value - mod.value)/err.value)

    @staticmethod
    def deriv_func(params, *args, **kwargs):
        """Jacobian of *fit_func*."""
        model = args[0]
        x = kwargs['x']
        y = kwargs['y']
        z = kwargs.get('z')
        if z is not None:
            err = kwargs.get('err', 1.*z.unit).to(z.unit)
            derivs = np.array([der.value for der in model.deriv(x, y, z.unit,
                *params)])
        else:
            err = kwargs.get('err', 1.*y.unit).to(y.unit)
            derivs = np.array([der.value for der in model.deriv(x, y.unit,
                *params)])

        diff = [np.ravel(chi) for chi in derivs/np.ravel(err.value)]
        return np.array(diff).T

