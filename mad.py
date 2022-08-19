import numpy as np


def _contains_nan(a, nan_policy='propagate'):
    policies = ['propagate', 'raise', 'omit']
    if nan_policy not in policies:
        raise ValueError("nan_policy must be one of {%s}" %
                         ', '.join("'%s'" % s for s in policies))
    try:
        # Calling np.sum to avoid creating a huge array into memory
        # e.g. np.isnan(a).any()
        with np.errstate(invalid='ignore'):
            contains_nan = np.isnan(np.sum(a))
    except TypeError:
        # This can happen when attempting to sum things which are not
        # numbers (e.g. as in the function `mode`). Try an alternative method:
        try:
            contains_nan = np.nan in set(a.ravel())
        except TypeError:
            # Don't know what to do. Fall back to omitting nan values and
            # issue a warning.
            contains_nan = False
            nan_policy = 'omit'
            warnings.warn("The input array could not be properly checked for nan "
                          "values. nan values will be ignored.", RuntimeWarning)

    if contains_nan and nan_policy == 'raise':
        raise ValueError("The input contains nan values")

    return (contains_nan, nan_policy)






def median_absolute_deviation(x, axis=0, center=np.median, scale=1.4826,
                              nan_policy='propagate'):
    """
    Compute the median absolute deviation of the data along the given axis.
    The median absolute deviation (MAD, [1]_) computes the median over the
    absolute deviations from the median. It is a measure of dispersion
    similar to the standard deviation, but is more robust to outliers [2]_.
    The MAD of an empty array is ``np.nan``.
    .. versionadded:: 1.3.0
    Parameters
    ----------
    x : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the range is computed. Default is 0. If None, compute
        the MAD over the entire array.
    center : callable, optional
        A function that will return the central value. The default is to use
        np.median. Any user defined function used will need to have the function
        signature ``func(arr, axis)``.
    scale : int, optional
        The scaling factor applied to the MAD. The default scale (1.4826)
        ensures consistency with the standard deviation for normally distributed
        data.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate'
        returns nan, 'raise' throws an error, 'omit' performs the
        calculations ignoring nan values. Default is 'propagate'.
    Returns
    -------
    mad : scalar or ndarray
        If ``axis=None``, a scalar is returned. If the input contains
        integers or floats of smaller precision than ``np.float64``, then the
        output data-type is ``np.float64``. Otherwise, the output data-type is
        the same as that of the input.
    See Also
    --------
    numpy.std, numpy.var, numpy.median, scipy.stats.iqr, scipy.stats.tmean,
    scipy.stats.tstd, scipy.stats.tvar
    Notes
    -----
    The `center` argument only affects the calculation of the central value
    around which the MAD is calculated. That is, passing in ``center=np.mean``
    will calculate the MAD around the mean - it will not calculate the *mean*
    absolute deviation.
    References
    ----------
    .. [1] "Median absolute deviation" https://en.wikipedia.org/wiki/Median_absolute_deviation
    .. [2] "Robust measures of scale" https://en.wikipedia.org/wiki/Robust_measures_of_scale
    Examples
    --------
    When comparing the behavior of `median_absolute_deviation` with ``np.std``,
    the latter is affected when we change a single value of an array to have an
    outlier value while the MAD hardly changes:
    >>> from scipy import stats
    >>> x = stats.norm.rvs(size=100, scale=1, random_state=123456)
    >>> x.std()
    0.9973906394005013
    >>> stats.median_absolute_deviation(x)
    1.2280762773108278
    >>> x[0] = 345.6
    >>> x.std()
    34.42304872314415
    >>> stats.median_absolute_deviation(x)
    1.2340335571164334
    Axis handling example:
    >>> x = np.array([[10, 7, 4], [3, 2, 1]])
    >>> x
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> stats.median_absolute_deviation(x)
    array([5.1891, 3.7065, 2.2239])
    >>> stats.median_absolute_deviation(x, axis=None)
    2.9652
    """
    x = np.asarray(x)

    # Consistent with `np.var` and `np.std`.
    if not x.size:
        return np.nan

    contains_nan, nan_policy = _contains_nan(x, nan_policy)

    if contains_nan and nan_policy == 'propagate':
        return np.nan

    if contains_nan and nan_policy == 'omit':
        # Way faster than carrying the masks around
        arr = ma.masked_invalid(x).compressed()
    else:
        arr = x

    if axis is None:
        med = center(arr)
        mad = np.median(np.abs(arr - med))
    else:
        med = np.apply_over_axes(center, arr, axis)
        mad = np.median(np.abs(arr - med), axis=axis)

    return scale * mad
