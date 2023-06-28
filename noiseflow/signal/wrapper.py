import os
import json
import numpy as np

from noiseflow.signal.python.taper import taper_py
from noiseflow.signal.python.detrend import detrend_py
from noiseflow.signal.python.decimate import decimate_py
from noiseflow.signal.python.filter import bandpass_py, bandstop_py, lowpass_py, highpass_py

# Determine the absolute path to the config file
env = os.environ.get('CONDA_DEFAULT_ENV')
config_path = os.path.abspath(os.path.expanduser(f'~/.noiseflow_config_{env}.json'))    

# Read the config.json file
with open(config_path) as f:
    config = json.load(f)

if config.get('NOISEFLOW_USE_CPP', False):
    NOISEFLOW_USE_CPP = True
    from noiseflow.lib import signal_share
else:
    NOISEFLOW_USE_CPP = False 



def detrend(data, type='linear', flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    # check type
    if type not in ['linear', 'demean']:
        raise ValueError("type must be 'linear' or 'demean'")
    
    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            detrend_py(data, type, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.detrend_float(data, type, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.detrend_double(data, type, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    

    
def taper(data, max_percentage=0.05, type='hann', side='both', flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check max_percentage
    if max_percentage < 0 or max_percentage > 0.5:
        raise ValueError("max_percentage must be between 0 and 0.5")
    
    # check type
    if type not in ['hann', 'hamming', 'blackman', 'blackman_harris', 'gaussian', 'triangular', 'bartlett', 'cosine', 'cosine_np', 'bartlett_hann', 'bohman', 'lanczos', 'flattop', 'kaiser']:
        raise ValueError("type must be 'hann', 'hamming', 'blackman', 'blackman_harris', 'gaussian', 'triangular', 'bartlett', 'cosine', 'cosine_np', 'bartlett_hann', 'bohman', 'lanczos', 'flattop', 'kaiser'")

    # check side
    if side not in ['both', 'left', 'right']:
        raise ValueError("side must be 'both', 'left', or 'right'")
    
    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            taper_py(data, max_percentage, type, side, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.taper_float(data, max_percentage, type, side, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.taper_double(data, max_percentage, type, side, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    


def bandpass(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check freq and df
    if (freqmin >= df/2) or (freqmax >= df/2):
        raise ValueError("freq must be < df/2")

    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            bandpass_py(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.bandpass_float(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.bandpass_double(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    


def bandstop(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check freq and df
    if (freqmin >= df/2) or (freqmax >= df/2):
        raise ValueError("freq must be < df/2")

    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            bandstop_py(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.bandstop_float(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.bandstop_double(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    


def lowpass(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check freq and df
    if freq >= df / 2:
        raise ValueError("freq must be < df/2")

    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            lowpass_py(data, freq, df, corners, zerophase, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.lowpass_float(data, freq, df, corners, zerophase, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.lowpass_double(data, freq, df, corners, zerophase, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    


def highpass(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check freq and df
    if freq >= df / 2:
        raise ValueError("freq must be < df/2")

    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])
    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            highpass_py(data, freq, df, corners, zerophase, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            signal_share.highpass_float(data, freq, df, corners, zerophase, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            signal_share.highpass_double(data, freq, df, corners, zerophase, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    


def decimate(data, df, df_new, flag=False, flag_gap=None, threads=1, py=False):
    # check demision of data
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # check df and df_new
    if df % df_new != 0:
        raise ValueError("Note df % df_new == 0")
    if df_new > df:
        raise ValueError("df_new must be <= df")
    
    # check threads
    if threads > data.shape[0]:
        raise ValueError("threads must be <= (data.shape[0]=%d)" % data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(data.shape[0])

    if py:
        if data.dtype == np.dtype(np.float32) or data.dtype == np.dtype(np.float64):
            results = decimate_py(data, df, df_new, flag, flag_gap)
        else:
            raise TypeError("data.dtype must be np.float32 or np.float64")
    else:
        if data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            results = signal_share.decimate_float(data, df, df_new, flag, flag_gap, threads)
        elif data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            results = signal_share.decimate_double(data, df, df_new, flag, flag_gap, threads)
        else:
            raise TypeError("data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    
    return results

