import scipy
import numpy as np


def decimate_py(data, df, df_new, flag=False, flag_gap=None):
    if data.ndim == 1:
        data = data.reshape(1, -1)

    test = decimate(data[0,:], df, df_new)
    redata = np.empty((data.shape[0], test.shape[0]))

    for i in range(0, data.shape[0]):
        redata[i,:] = decimate(data[i,:], df, df_new)

    return redata


def decimate(data, df, df_new):
    factor = int(df / df_new)
    if factor > 13:
        msg = (
                "IRR filter is unstable for decimation factors above"
                " 13. Call decimate multiple times."
            )
        raise ValueError(msg)
    
    redata = scipy.signal.decimate(data, factor, ftype='iir', axis=-1, zero_phase=True)

    return redata




