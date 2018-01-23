import numpy as np
from scipy.stats import binned_statistic

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
            it[0] = functions[kwargs['statistic']](values.flatten()[ind])
        except ZeroDivisionError:
            it[0] = np.nan
        it.iternext()

    return out

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

