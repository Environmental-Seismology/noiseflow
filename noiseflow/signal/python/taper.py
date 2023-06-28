import scipy
import warnings
import numpy as np


def taper_py(data, max_percentage=0.05, type='hann', side='both', flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    for i in range(0, data.shape[0]):
        taper(data[i,:], max_percentage, type, side)



def taper(data, max_percentage=0.05, type='hann', side='both'):

    side_valid = ['both', 'left', 'right']
    npts = len(data)
    if side not in side_valid:
        raise ValueError("'side' has to be one of: %s" % side_valid)

    # max_half_lenghts_1
    max_half_lenghts_1 = None
    if max_percentage is not None:
        max_half_lenghts_1 = int(max_percentage * npts)

    if 2 * max_half_lenghts_1 > npts:
        msg = "The requested taper is longer than the data. " \
                "The taper will be shortened to data length."
        warnings.warn(msg)

    # max_half_lenghts_2
    max_half_lenghts_2 = int(npts / 2)
    if max_half_lenghts_1 is None:
        wlen = max_half_lenghts_2
    else:
        wlen = min(max_half_lenghts_1, max_half_lenghts_2)

    if type == "hann":
        taper_sides = scipy.signal.windows.hann(2*wlen+1)
    elif type == "hamming":
        taper_sides = scipy.signal.windows.hamming(2*wlen+1)
    elif type == "bartlett":
        taper_sides = scipy.signal.windows.bartlett(2*wlen+1)
    elif type == "blackman":
        taper_sides = scipy.signal.windows.blackman(2*wlen+1)
    elif type == "flattop":
        taper_sides = scipy.signal.windows.flattop(2*wlen+1)
    elif type == "parzen":
        taper_sides = scipy.signal.windows.parzen(2*wlen+1)
    elif type == "bohman":
        taper_sides = scipy.signal.windows.bohman(2*wlen+1)
    elif type == "blackmanharris":
        taper_sides = scipy.signal.windows.blackmanharris(2*wlen+1)
    elif type == "nuttall":
        taper_sides = scipy.signal.windows.nuttall(2*wlen+1)
    elif type == "barthann":
        taper_sides = scipy.signal.windows.barthann(2*wlen+1)
    elif type == "kaiser":
        taper_sides = scipy.signal.windows.kaiser(2*wlen+1, beta=14)
    elif type == "gaussian":
        taper_sides = scipy.signal.windows.gaussian(2*wlen+1, std=14)
    elif type == "general_gaussian":
        taper_sides = scipy.signal.windows.general_gaussian(2*wlen+1, p=1, sig=14)
    else:
        raise ValueError("Unknown taper type: %s" % type)

    if side == 'left':
        taper = np.hstack((taper_sides[:wlen], np.ones(npts - wlen)))
    elif side == 'right':
        taper = np.hstack((np.ones(npts - wlen),
                            taper_sides[len(taper_sides) - wlen:]))
    else:
        taper = np.hstack((taper_sides[:wlen], np.ones(npts - 2 * wlen),
                            taper_sides[len(taper_sides) - wlen:]))

    # Convert data if it's not a floating point type.
    if not np.issubdtype(data.dtype, np.floating):
        data = np.require(data, dtype=np.float64)

    data *= taper




