import numpy as np
from astropy.modeling import fitting

def bootstrap(x, y, err, p_init, niter=100, fit='Guassian'):
    res = []
    for i in range(niter):
        # Generate random points
        rand_data = np.array([np.random.normal(yi, erri, 1)[0] if erri!=0 else yi for yi, erri in zip(y,err)])

        # Fit
        if fit.lower()=='gaussian':
            fit_p = fitting.LevMarLSQFitter()
            p = fit_p(p_init, x, rand_data)
            res += [p.stddev.value]
        elif fit.lower()=='linear':
            fit_p = fitting.LinearLSQFitter()
            p = fit_p(p_init, x, rand_data)
            res += [[p.slope.value, p.intercept.value]]
        elif fit.lower()=='power':
            fit_p = fitting.LevMarLSQFitter()
            p = fit_p(p_init, x, rand_data)
            res += [[p.amplitude.value, p.alpha.value, p.x_0.value]]
        else:
            raise ValueError('Fit (%s) not implemented' % fit)
                    
    return np.array(res)

def resample(data, error, distribution='normal'):
    vfunc = np.vectorize(np.random.normal)
    return vfunc(data, error)
