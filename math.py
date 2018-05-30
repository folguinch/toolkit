import numpy as np
from scipy.stats import binned_statistic
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d

def rms(x):
    """Root mean square (rms).

    The rms of *x* is defined as:
    rms = sqrt( sum(x_i**2) / n )
    where n is the number of points in x.

    Parameters:
        x (np.array): data.

    Return:
        rms (float): the root mean square of *x*.
    """
    assert hasattr(x, 'size')
    return np.sqrt(np.sum(x**2)/x.size)

def chi2(obs, mod, err=None, mask=None, dof=1, logger=None):
    """Calculate the chi2 between two arrays.

    If errors are not given, the chi2 is equivalent to the mean average
    difference (mad) squared when dof=0.

    The mad is defined as:
        mad = sum(|obs_i - mod_i|)/N
    with N the number of valid points. The chi2 is defined as:
        chi2 = sum((obs_i - mod_i)**2 / err_i**2)/(N - dof)

    The mask should be defined so it is True where the data is compared. The
    number of valid points N is equal to the sum of the True values in the mask.

    Parameters:
        obs (np.array): observed value.
        mod (np.array): expected value.
        err (np.array, default=None): errors.
        mask (np.array, default=None): mask for the data.
        dof (int, default=0): generally, the number of parameters to fit.

    Return:
        mad, chi (floats): the mad and chi2 values.
    """

    # Create a default mask
    if mask is None:
        mask = np.ones(obs.shape, dtype=bool)
    assert obs.shape==mod.shape==mask.shape

    # Create error array 
    if err is None:
        err = np.ones(obs.shape)
    elif not hasattr(err, 'shape'):
        err = np.ones(obs.shape) * err

    # Statistics
    mad = np.abs(obs - mod)
    chi2 = ((obs-mod)/err)**2

    # Number of points
    N = np.sum(mask)
    if logger:
        logger.info('Valid points: %i / %i', N, chi2.size)
    else:
        print 'Valid points: %i / %i' % (N, chi2.size)

    return np.mean(mad[mask]), np.sum(chi2[mask])/float(N - dof)

#def round_to_sigfig(x, n):
#    return np.around(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def binned_stat(values, x, statistic='mean', bins=10, range=None):
    """Redefine np.binned_statistic for np.apply_along_axis.

    Parameters:
        Same as np.binned_statistic but invert x and value positions.

    Returns:
        The binned values only.        
    """
    return binned_statistic(x, values, statistic=statistic, bins=bins,
            range=range)[0]

def rebin_along_axis(x, values, axis, statistic='sum', bins=10, range=None):
    return np.apply_along_axis(binned_stat, axis, values, x, 
            statistic=statistic, bins=bins, range=range)

def rebin_regular_nd(values, *args, **kwargs):
    """Rebin a n-dimensional array.

    The input array should be in a regular array, e.g. images. 

    Parameters:
        values: the data values.
        args: the coordinates along each axis to be binned.
        kwargs: the available keys are:
            * bins: the number of bins or the bin edges.
            * statistic (default=mean): the statistic to apply (mean, sum,
                count).
    """
    # Checks
    av_keys = ['bins', 'statistic']
    assert all(map(lambda x: x in av_keys, kwargs))
    assert kwargs.setdefault('statistic', 'mean') in ['mean', 'sum', 'count']

    sums = values
    counts = np.ones(values.shape)
    for i in range(values.ndim):
        try:
            bins = kwargs.get('bins', [10]*values.ndim)[i]
        except TypeError:
            bins = kwargs.get('bins')
        if kwargs['statistic'] in ['mean', 'sum']:
            sums = rebin_along_axis(args[i], sums, i, bins=bins)
        if kwargs['statistic'] in ['mean', 'count']:
            counts = rebin_along_axis(args[i], counts, i, bins=bins)

    if kwargs['statistic']=='mean':
        return sums/counts
    elif kwargs['statistic']=='sum':
        return sums
    else:
        return counts

def rebin_irregular_nd(values, bins, *args, **kwargs):
    """Rebin a n-dimesnion array.

    The input array can be defined in an irregular array, e.g. by converting
    regular spherical to irregular cartesian coordinates.

    The output array will have a shape *len(bin)-1* for each bin in bins.

    Parameters:
        values: the data values.
        bins: the bin edges.
        args: coordinates along each axis to be binned. Each corrdinate must
            have the same shape as values.
        kwargs: available keys are:
            * statistic (default=mean): statistic to calculate over the data
                values (mean, average, sum, count).
            * weights: weights for the average statistic.
    """
    # Checks
    functions = {'mean': np.mean, 'average': np.average, 'sum': np.sum} #, 'count': count}
    av_keys = ['statistic', 'weights']
    assert all(map(lambda x: x in av_keys, kwargs))
    assert kwargs.setdefault('statistic', 'mean') in functions
    assert len(bins)==values.ndim

    # Digitize the axes
    bins_d = []
    new_shape = []
    for i in range(values.ndim):
        # Emulate the behaviour of scipy binned_statistic, i.e. include the
        # last bin edge in the binning.
        dig = np.digitize(args[i].flatten(), bins[i])
        dig[args[i].flatten() == bins[i][-1]] = len(bins[i])-1

        # Append to list
        bins_d += [dig]
        new_shape += [len(bins[i])-1]

    # Fill the new grid
    out = np.zeros(new_shape)
    it = np.nditer(out, flags=['multi_index'], op_flags=['readwrite'])
    while not it.finished:
        ind = np.ones(len(bins_d[0]), dtype=bool)
        for i,b in enumerate(bins_d):
            ind = np.logical_and(ind, b==it.multi_index[i]+1)

        try:
            it[0] = functions[kwargs['statistic']](values.flatten()[ind],
                    weights=kwargs.get('weights').flatten()[ind])
        except TypeError:
            print 'Function does not have weight keyword'
            it[0] = functions[kwargs['statistic']](values.flatten()[ind])
        except ZeroDivisionError:
            it[0] = np.nan
        it.iternext()

    return out

def map_sph_to_cart(val, new_x, new_y, new_z, r, th, phi=None, **kwargs):
    """Resample the distribution in val(r, th, phi) in new Cartesian
    coordinates.

    Solution taken from:
    https://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid

    Parameters:
        val: the values at the input position.
        new_x, new_y, new_z: coordinates of the new grid points.
        r, th: coordinates of the values in spherical coordinates.
        phi (default=None): azimuthal spherical coordinate angle.
        kwargs: keyword arguments for scipy map_coordinates function.

    Notes:
        The default `order` keyword for scipy map_coordinates is set to zero,
        i.e. nearest neighbour.
    """
    # Set default interpolation order to nearest neighbour
    kwargs.setdefault('order', 0)

    # New grid meshes
    if new_x.ndim==3 and new_y.ndim==3 and new_z.ndim==3:
        Z, Y, X = new_z, new_y, new_x
    else:
        Z, Y, X = np.meshgrid(new_z, new_y, new_x, indexing='ij')
    X = X.to(r.unit).value
    Y = Y.to(r.unit).value
    Z = Z.to(r.unit).value

    # New grid in spherical coordinates
    new_r = np.sqrt(X**2 + Y**2 + Z**2)
    new_th = np.arctan2(np.sqrt(X**2+Y**2), Z)
    
    # Coordinate interpolation functions
    ir = interp1d(r.value, np.arange(len(r)), bounds_error=False)
    ith = interp1d(th.value, np.arange(len(th)), bounds_error=False)

    # New coordinates indices
    new_ir = ir(new_r.ravel())
    new_ith = ith(new_th.ravel())

    # Values outside the grid
    new_ir[new_r.ravel() > r.value.max()] = len(r)-1
    new_ir[new_r.ravel() < r.value.min()] = 0
    new_ith[new_th.ravel() > th.value.max()] = len(th)-1
    new_ith[new_th.ravel() < th.value.min()] = 0

    # Include phi if needed
    if phi is not None:
        new_phi = np.arctan2(Y, X)
        iphi = interp1d(phi.value, np.arange(len(phi)), bounds_error=False)
        new_iphi = iphi(new_phi.ravel())
        new_iphi[new_phi.ravel() > phi.value.max()] = len(phi)-1
        new_iphi[new_phi.ravel() < phi.value.min()] = 0
        map_to = [new_iphi, new_it, new_ir]
    else:
        map_to = [new_ith, new_ir]

    # Output
    out = map_coordinates(val, map_to, **kwargs)

    # Values outside the largest radius set to zero
    out[new_r.ravel() > r.value.max()] = 0.

    return out.reshape(new_r.shape)

def map_sph_to_cart_axisym(val, r, th, new_x, new_y, new_z, **kwargs):
    """Resample the axisymmetric distribution in val(r, th) into new Cartesian
    coordinates.

    Parameters:
        val: the values at the input positions.
        r, th: 1-D array of the *val* coordinates.
        new_x, new_y, new_z: Cartesian coordinates of the new grid points.
        kwargs: keyword arguments for scipy map_coordinates function.

    Notes:
        The default `order` keyword for scipy map_coordinates is set to zero,
        i.e. nearest neighbour.
    """
    assert val.ndim == 2
    return map_sph_to_cart(val, new_x, new_y, new_z, r, th, **kwargs)

def map_sph_to_cart_3d(val, r, th, phi, new_x, new_y, new_z, **kwargs):
    """Resample the distribution in val(r, th, phi) into new Cartesian
    coordinates.

    This is the same as the `map_sph_to_cart` function above but with a
    different input order and checks.

    Parameters:
        val: the values at the input positions.
        r, th, phi: 1-D array of the *val* coordinates.
        new_x, new_y, new_z: Cartesian coordinates of the new grid points.
        kwargs: keyword arguments for scipy map_coordinates function.

    Notes:
        The default `order` keyword for scipy map_coordinates is set to zero,
        i.e. nearest neighbour.
    """
    assert val.ndim == 3
    return map_sph_to_cart(val, new_x, new_y, new_z, r, th, phi=phi, **kwargs)

if __name__=='__main__':
    z = np.array([[1,2,3,4,5,6],
        [7,8,9,10,11,12],
        [12,11,10,9,8,7],
        [6,5,4,3,2,1],
        [1,2,3,4,5,6],
        [1,2,3,4,5,6]])

    x = np.arange(6)
    X, Y = np.meshgrid(x,x)
    bins = [np.array([ 0.        ,  1.66666667,  3.33333333,  5.        ]),
            np.array([ 0.        ,  1.66666667,  3.33333333,  5.        ])]

    print z

    print rebin_regular_nd(z, x, x, statistic='mean', bins=3)
    print rebin_irregular_nd(z, bins, Y, X, statistic='mean')

