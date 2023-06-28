import numpy as np

from scipy.signal import detrend as scipy_detrend
from scipy.interpolate import LSQUnivariateSpline


def detrend_py(data, type='linear', flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    for i in range(0, data.shape[0]):
        detrend(data[i,:], type)


def spline(data, order=2, dspline=1000):
    if not np.issubdtype(data.dtype, np.floating):
        data = np.require(data, dtype=np.float64)

    x = np.arange(len(data))
    splknots = np.arange(dspline/2.0, len(data) - dspline/2.0 + 2, dspline)
    spl = LSQUnivariateSpline(x=x, y=data, t=splknots, k=order)
    fit = spl(x)
    data -= fit

    return data


def detrend(data, type='linear', order=2, dspline=1000):
    """
    Detrend data.
    :param data: The data to detrend.
    :type data: :class:`numpy.ndarray`
    :param type: The type of detrending to perform. Can be one of ``'demean'``,
        ``'linear'``, or ``'spline'``.
    :type type: str
    :param overwrite_data: If True the data will be detrended in place.
    :type overwrite_data: bool
    :param dtype: The data type to use for the detrended data. Defaults to
        ``np.float32``.
    :type dtype: :class:`numpy.dtype`
    :param order: The order of the spline to fit. Defaults to 2.
    :type order: int
    :param dspline: The distance between spline knots. Defaults to 1000.
    :type dspline: int

    :returns: The detrended data.
    :rtype: :class:`numpy.ndarray`

    .. rubric:: Example

    >>> import numpy as np
    >>> from noiseflow.signal import detrend
    >>> data = np.array([1.1,2.7,3,3,3,6,7,8,9,10], overwrite_data=True)
    >>> detrend(data, type='spline')

    """

    if type == 'demean':
        data -= np.mean(data)
    elif type == 'linear':
        scipy_detrend(data, type='linear', overwrite_data=True)
    elif type == 'spline':
        spline(data, order=order, dspline=dspline)
    else:
        raise ValueError("type must be 'demean', 'linear', or 'spline'")






