import warnings

from scipy.signal import iirfilter
from scipy.signal import sosfilt
from scipy.signal import zpk2sos



def bandpass_py(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    for i in range(0, data.shape[0]):
        data[i,:] = bandpass(data[i,:], freqmin, freqmax, df, corners, zerophase)


def bandstop_py(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    for i in range(0, data.shape[0]):
        data[i,:] = bandstop(data[i,:], freqmin, freqmax, df, corners, zerophase)


def lowpass_py(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    for i in range(0, data.shape[0]):
        data[i,:] = lowpass(data[i,:], freq, df, corners, zerophase)


def highpass_py(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)
        
    for i in range(0, data.shape[0]):
        data[i,:] = highpass(data[i,:], freq, df, corners, zerophase)



def bandpass(data, freqmin, freqmax, df, corners=4, zerophase=False):
    """
    Butterworth-Bandpass Filter.
    Filter data from ``freqmin`` to ``freqmax`` using ``corners``
    corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).
    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freqmin: Pass band low corner frequency.
    :param freqmax: Pass band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the filter order but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high - 1.0 > -1e-6:
        msg = ("Selected high corner frequency ({}) of bandpass is at or "
               "above Nyquist ({}). Applying a high-pass instead.").format(
            freqmax, fe)
        warnings.warn(msg)
        return highpass(data, freq=freqmin, df=df, corners=corners,
                        zerophase=zerophase)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, [low, high], btype='band',
                        ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)


def bandstop(data, freqmin, freqmax, df, corners=4, zerophase=False):
    """
    Butterworth-Bandstop Filter.
    Filter data removing data between frequencies ``freqmin`` and ``freqmax``
    using ``corners`` corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).
    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freqmin: Stop band low corner frequency.
    :param freqmax: Stop band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, [low, high],
                        btype='bandstop', ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)


def lowpass(data, freq, df, corners=4, zerophase=False):
    """
    Butterworth-Lowpass Filter.
    Filter data removing data over certain frequency ``freq`` using ``corners``
    corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).
    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freq: Filter corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    f = freq / fe
    # raise for some bad scenarios
    if f > 1:
        f = 1.0
        msg = "Selected corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    z, p, k = iirfilter(corners, f, btype='lowpass', ftype='butter',
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)


def highpass(data, freq, df, corners=4, zerophase=False):
    """
    Butterworth-Highpass Filter.
    Filter data removing data below certain frequency ``freq`` using
    ``corners`` corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).
    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freq: Filter corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    f = freq / fe
    # raise for some bad scenarios
    if f > 1:
        msg = "Selected corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, f, btype='highpass', ftype='butter',
                        output='zpk')
    sos = zpk2sos(z, p, k)
    if zerophase:
        firstpass = sosfilt(sos, data)
        return sosfilt(sos, firstpass[::-1])[::-1]
    else:
        return sosfilt(sos, data)






