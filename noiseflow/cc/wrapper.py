import os
import re
import json
import numpy as np

from noiseflow.cc.rfftdata import RFFTData_Class
from noiseflow.cc.corrdata import CorrData_Class
from noiseflow.cc.stackdata import StackData_Class
from noiseflow.cc.python.rfft import RFFTClass_python
from noiseflow.cc.python.corr import CorrClass_python
from noiseflow.cc.python.stack import StackClass_python

# Determine the absolute path to the config file
env = os.environ.get('CONDA_DEFAULT_ENV')
env = re.sub('[^a-zA-Z0-9_]', '', env)[0:50]
config_path = os.path.abspath(os.path.expanduser(f'~/.noiseflow_config_{env}.json'))    

# Read the config.json file
with open(config_path) as f:
    config = json.load(f)

if config.get('NOISEFLOW_USE_CPP', False):
    NOISEFLOW_USE_CPP = True
    from noiseflow.lib import cc_share
else:
    NOISEFLOW_USE_CPP = False  


def rfft(raw_data, dt,
        cc_len, cc_step, 
        time_norm, clip_std, smooth_N,
        freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N, 
        flag=False, flag_gap=None, threads=1, jobs=1, py=False):
    
    # check demision of raw_data
    if raw_data.ndim == 1:
        raw_data = raw_data.reshape(1, -1)

    # check freqmax
    if (freqmax > (1/(2*dt))):
        raise ValueError("freqmax must be <= (%fhz) according to fmax=1/(2*dt) formula" % (1/(2*dt)))
    
    # check cc_len
    if (cc_len > raw_data.shape[1]*dt):
        raise ValueError("cc_len must be <= (raw_data.shape[1]*dt=%f)" % (raw_data.shape[1]*dt))

    # check cc_step
    if (cc_step >= cc_len):
        raise ValueError("cc_step must be < (cc_len=%f)" % cc_len)

    # check time_norm
    if time_norm not in ['no', 'onebit', 'clip', 'smooth']:
        raise ValueError("time_norm must be 'no', 'onebit', 'clip', or 'smooth'")
    
    # check freq_norm
    if freq_norm not in ['no', 'whiten', 'smooth_whiten']:
        raise ValueError("freq_norm must be 'no', 'whiten', or 'smooth_whiten'")
    
    # check threads
    if threads > raw_data.shape[0]:
        raise ValueError("threads must be <= (raw_data.shape[0]=%d)" % raw_data.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(raw_data.shape[0])

    # run rfft
    if py:
        if raw_data.dtype == np.dtype(np.float32) or raw_data.dtype == np.dtype(np.float64):
            m = RFFTClass_python(raw_data, dt,
                         cc_len, cc_step, 
                        time_norm, clip_std, smooth_N,
                        freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
                        flag, jobs)
            m.run()
            rfft_data = m.output_data
        else:
            raise TypeError("raw_data.dtype must be np.float32 or np.float64")
    else:
        if raw_data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            rfft_data = cc_share.rfft_float(raw_data, dt, 
                cc_len, cc_step, 
                time_norm, clip_std, smooth_N,
                freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
                flag, flag_gap,
                threads)
        elif raw_data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            rfft_data = cc_share.rfft_double(raw_data, dt, 
                cc_len, cc_step, 
                time_norm, clip_std, smooth_N,
                freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
                flag, flag_gap,
                threads)
        else:
            raise TypeError("raw_data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    
    # generate RFFTData object
    RFFTData = RFFTData_Class(rfft_data, dt, cc_len, cc_step, 
        time_norm, clip_std, smooth_N,
        freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
        flag, flag_gap,
        threads, jobs, py)
    
    return RFFTData



def corr(rfft_data, dt, corr_method, corr_pair, maxlag, smoothspect_N=10, 
         flag=False, flag_gap=None, threads=1, jobs=1, py=False):
    
    # check corr_method
    if corr_method not in ['xcorr', 'deconv', 'coherency']:
        raise ValueError("corr_method must be 'xcorr', 'deconv', or 'coherency'")
    
    # check maxlag
    if maxlag > dt*rfft_data.shape[2]:
        raise ValueError("maxlag must be <= (dt*rfft_data.shape[2]=%f)" % (dt*rfft_data.shape[2]))
    
    # check threads
    if threads > rfft_data.shape[0]:
        raise ValueError("threads must be <= (rfft_data.shape[0]=%d)" % rfft_data.shape[0])

    # check demision of corr_pair
    if corr_pair.ndim == 1:
        corr_pair = corr_pair.reshape(1, -1)

    # check threads
    if threads > corr_pair.shape[0]:
        raise ValueError("threads must be <= (corr_pair.shape[0]=%d)" % corr_pair.shape[0])
    
    # check flag_gap
    if flag_gap==None:
        flag_gap = int(corr_pair.shape[0])

    # run corr
    if py:
        if rfft_data.dtype == np.dtype(np.complex64) or rfft_data.dtype == np.dtype(np.complex128):
            m = CorrClass_python(rfft_data, dt, corr_method, corr_pair, maxlag,
                smoothspect_N, flag, jobs)
            m.run()
            corr_data = m.output_data
        else:
            raise TypeError("rfft_data.dtype must be np.complex64 or np.complex128")
    else:
        if rfft_data.dtype == np.dtype(np.complex64) and NOISEFLOW_USE_CPP:
            corr_data = cc_share.corr_float(rfft_data, dt, corr_method, corr_pair, maxlag,
                smoothspect_N, flag, flag_gap, threads)
        elif rfft_data.dtype == np.dtype(np.complex128) and NOISEFLOW_USE_CPP:
            corr_data = cc_share.corr_double(rfft_data, dt, corr_method, corr_pair, maxlag,
                smoothspect_N, flag, flag_gap, threads)
        else:
            raise TypeError("rfft_data.dtype must be np.complex64, np.complex128 or NOISEFLOW_USE_CPP=False")

    # generate CorrData object
    CorrData = CorrData_Class(corr_data, dt, corr_method, corr_pair, maxlag, smoothspect_N,
        flag, flag_gap, threads, jobs, py)
    
    return CorrData



def stack(corr_data, dt, stack_method, 
          par=None, 
          stack_all=True, stack_len=0, stack_step=0, 
          pick=False, median_high=10, median_low=0.1, 
          flag=False, flag_gap=None, threads=1, jobs=1, py=False):
    
    # check par
    if par == None:
        par = {"axis":0,"p":2,"g":1,"cc_min":0.0,"epsilon":1E-5,"maxstep":10,
                "win":None,"stat":False,"h":0.75,'plot':False,'normalize':True,'ref':None},
    
    # check demision of corr_data
    if corr_data.ndim == 2:
        corr_data = corr_data.reshape(1, corr_data.shape[0], corr_data.shape[1])

    # check stack_method
    if stack_method not in ["linear","pws","robust","acf","nroot","selective","cluster","tfpws"]:
        raise ValueError("stack_method must be 'linear', 'pws', 'robust', 'acf', 'nroot', 'selective', 'cluster', or 'tfpws'")
    
    # check stack_len
    if (stack_len > corr_data.shape[1]):
        raise ValueError("stack_len must be <= (corr_data.shape[1]=%d)" % corr_data.shape[1])

    # check stack_step
    if (stack_step >= stack_len):
        raise ValueError("stack_step must be < (stack_len=%d)" % stack_len)
    
    # check threads
    if threads > corr_data.shape[0]:
        raise ValueError("threads must be <= (corr_data.shape[0]=%d)" % corr_data.shape[0])

    # check flag_gap
    if flag_gap==None:
        flag_gap = int(corr_data.shape[0])

    # run stack
    if py:
        if corr_data.dtype == np.dtype(np.float32) or corr_data.dtype == np.dtype(np.float64):
            m = StackClass_python(corr_data, stack_method, par, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, jobs)
            m.run()
            stack_data = m.output_data
            stack_ngood = m.ngood_all
        else:
            raise TypeError("corr_data.dtype must be np.float32 or np.float64")
    else:
        if corr_data.dtype == np.dtype(np.float32) and NOISEFLOW_USE_CPP:
            m = cc_share.StackClass_float(corr_data, stack_method, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, flag_gap, threads)
            m.run()
            stack_data = m.get_stackdata()
            stack_ngood = m.get_ngood()
        elif corr_data.dtype == np.dtype(np.float64) and NOISEFLOW_USE_CPP:
            m = cc_share.StackClass_double(corr_data, stack_method, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, flag_gap, threads)
            m.run()
            stack_data = m.get_stackdata()
            stack_ngood = m.get_ngood()
        else:
            raise TypeError("corr_data.dtype must be np.float32, np.float64 or NOISEFLOW_USE_CPP=False")
    
    # generate StackData object
    StackData = StackData_Class(stack_data, stack_ngood, dt,
                 stack_method, par, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, flag_gap, threads, jobs, py)
    
    return StackData



