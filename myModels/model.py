import numpy as np
from astropy.modeling.fitting import LevMarLSQFitter, _convert_input, _model_to_fit_params

class Model(object):
    """Class to manage model information and parameters.

    Attributes:
        params (dict): model parameters with units.
    """

    def __init__(self, **params):
        self.params = params

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    def values(self, units=None):
        vals = {}
        for key,val in self.params.items():
            if units is not None and units.get(key):
                vals[key] = val.to(units[key]).value
            else:
                vals[key] = val
        return vals

    def fit(self, x, y, model, *args):
        """Fit the model to the given data.

        Parameters:
            x: x-values.
            y: values to fit.
            model: astropy model function.
        """
        p_init = model(*args, **self.values(model.STD_UNITS))
        fit_p = myLevMarLSQFitter()
        p = fit_p(p_init, x.to(model.STD_UNITS['x']).value,
                y.to(model.STD_UNITS['y']).value)
        print p


class myLevMarLSQFitter(LevMarLSQFitter):
    def __call__(self, model, x, y, z=None, weights=None,
            maxiter=100, acc=1e-07,
            epsilon=np.sqrt(np.finfo(float).eps), 
            estimate_jacobian=False):
        from scipy import optimize
        farg = (model, weights, ) + _convert_input(x, y, z)
        print farg

        if model.fit_deriv is None or estimate_jacobian:
            dfunc = None
        else:
            dfunc = self._wrap_deriv
        init_values, _ = _model_to_fit_params(model)
        print init_values
        fitparams, cov_x, dinfo, mess, ierr = optimize.leastsq(
                self.objective_function, init_values, args=farg, Dfun=dfunc,
                col_deriv=model.col_fit_deriv, maxfev=maxiter, epsfcn=epsilon,
                xtol=acc, full_output=True)
        _fitter_to_model_params(model, fitparams)
        self.fit_info.update(dinfo)
        self.fit_info['cov_x'] = cov_x
        self.fit_info['message'] = mess
        self.fit_info['ierr'] = ierr
        if ierr not in [1, 2, 3, 4]:
            print "The fit may be unsuccessful; check fit_info['message'] for more information."
            #warnings.warn("The fit may be unsuccessful; check "
            #                "fit_info['message'] for more information.",
            #                AstropyUserWarning)

        # now try to compute the true covariance matrix
        if (len(y) > len(init_values)) and cov_x is not None:
            sum_sqrs = np.sum(self.objective_function(fitparams, *farg)**2)
            dof = len(y) - len(init_values)
            self.fit_info['param_cov'] = cov_x * sum_sqrs / dof
        else:
            self.fit_info['param_cov'] = None

        return model
