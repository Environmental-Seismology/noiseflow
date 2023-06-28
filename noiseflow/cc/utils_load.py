import h5py
import numpy as np

from scipy.io import loadmat
from noiseflow.signal.rawdata import RawData_Class
from noiseflow.cc.rfftdata import RFFTData_Class
from noiseflow.cc.corrdata import CorrData_Class
from noiseflow.cc.stackdata import StackData_Class


def load_raw(filename, format="npz", only_header=False):
    if format == "npz":
        DD = np.load(filename)
        if only_header:
            data = None
        else:
            data = DD["data"]

        sampling_rate = DD["sampling_rate"].reshape(1)[0]

    elif format == "h5":
        with h5py.File(filename, 'r') as f:
            if only_header:
                data = None
            else:
                data = f['noiseflow_group']['data'][:]

            sampling_rate = f['noiseflow_group'].attrs['sampling_rate']

    elif format == "mat":
        DD = loadmat(filename)
        if only_header:
            data = None
        else:
            data = DD["data"][:,:]

        sampling_rate = DD["sampling_rate"][0][0]

    else:
        ValueError("format must be 'npz', 'h5' or 'mat'")
    
    RawData = RawData_Class(data, sampling_rate)
    
    return RawData



def load_rfft(filename, format="npz", only_header=False):
    if format == "npz":
        DD = np.load(filename)
        if only_header:
            rfft_data = None
        else:
            rfft_data = DD["rfft_data"]

        dt = DD["dt"].reshape(1)[0]
        cc_len = DD["cc_len"].reshape(1)[0]
        cc_step = DD["cc_step"].reshape(1)[0]
        time_norm = DD["time_norm"].reshape(1)[0]
        clip_std = DD["clip_std"].reshape(1)[0]
        smooth_N = DD["smooth_N"].reshape(1)[0]
        freq_norm = DD["freq_norm"].reshape(1)[0]
        freqmin = DD["freqmin"].reshape(1)[0]
        freqmax = DD["freqmax"].reshape(1)[0]
        whiten_npad = DD["whiten_npad"].reshape(1)[0]
        smoothspect_N = DD["smoothspect_N"].reshape(1)[0]
        flag = DD["flag"].reshape(1)[0]
        flag_gap = DD["flag_gap"].reshape(1)[0]
        threads = DD["threads"].reshape(1)[0]
        jobs = DD["jobs"].reshape(1)[0]
        py = DD["py"].reshape(1)[0]

    elif format == "h5":
        with h5py.File(filename, 'r') as f:
            group = f['noiseflow_group']

            if only_header:
                rfft_data = None
            else:
                rfft_data = group['rfft_data'][:]

            dt = group.attrs['dt']
            cc_len = group.attrs['cc_len']
            cc_step = group.attrs['cc_step']
            time_norm = group.attrs['time_norm']
            clip_std = group.attrs['clip_std']
            smooth_N = group.attrs['smooth_N']
            freq_norm = group.attrs['freq_norm']
            freqmin = group.attrs['freqmin']
            freqmax = group.attrs['freqmax']
            whiten_npad = group.attrs['whiten_npad']
            smoothspect_N = group.attrs['smoothspect_N']
            flag = group.attrs['flag']
            flag_gap = group.attrs['flag_gap']
            threads = group.attrs['threads']
            jobs = group.attrs['jobs']
            py = group.attrs['py']

    elif format == "mat":
        DD = loadmat(filename)
        if only_header:
            rfft_data = None
        else:
            rfft_data = DD["rfft_data"][:,:,:]

        if DD["flag"][0][0] == 1:
            flag = True
        else:
            flag = False
        
        if DD["py"][0][0] == 1:
            py = True
        else:
            py = False

        dt = DD["dt"][0][0]
        cc_len = DD["cc_len"][0][0]
        cc_step = DD["cc_step"][0][0]
        time_norm = DD["time_norm"][0]
        clip_std = DD["clip_std"][0][0]
        smooth_N = DD["smooth_N"][0][0]
        freq_norm = DD["freq_norm"][0]
        freqmin = DD["freqmin"][0][0]
        freqmax = DD["freqmax"][0][0]
        whiten_npad = DD["whiten_npad"][0][0]
        smoothspect_N = DD["smoothspect_N"][0][0]
        flag_gap = DD["flag_gap"][0][0]
        threads = DD["threads"][0][0]
        jobs = DD["jobs"][0][0]

    else:
        ValueError("format must be 'npz', 'h5' or 'mat'")
    
    RFFTData = RFFTData_Class(rfft_data, dt, cc_len, cc_step, 
        time_norm, clip_std, smooth_N,
        freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
        flag, flag_gap,
        threads, jobs, py)
    
    return RFFTData



def load_corr(filename, format="npz", only_header=False):
    if format == "npz":
        DD = np.load(filename)
        if only_header:
            corr_data = None
        else:
            corr_data = DD["corr_data"]

        corr_pair = DD["corr_pair"]
        dt = DD["dt"].reshape(1)[0]
        corr_method = DD["corr_method"].reshape(1)[0]
        maxlag = DD["maxlag"].reshape(1)[0]
        smoothspect_N = DD["smoothspect_N"].reshape(1)[0]
        flag = DD["flag"].reshape(1)[0]
        flag_gap = DD["flag_gap"].reshape(1)[0]
        threads = DD["threads"].reshape(1)[0]
        jobs = DD["jobs"].reshape(1)[0]
        py = DD["py"].reshape(1)[0]

    elif format == "h5":
        with h5py.File(filename, 'r') as f:
            group = f['noiseflow_group']

            if only_header:
                corr_data = None
            else:
                corr_data = group['corr_data'][:]

            dt = group.attrs['dt']
            corr_method = group.attrs['corr_method']
            maxlag = group.attrs['maxlag']
            smoothspect_N = group.attrs['smoothspect_N']
            flag = group.attrs['flag']
            flag_gap = group.attrs['flag_gap']
            threads = group.attrs['threads']
            corr_pair = group['corr_pair'][:]
            jobs = group.attrs['jobs']
            py = group.attrs['py']

    elif format == "mat":
        DD = loadmat(filename)
        if only_header:
            corr_data = None
        else:
            corr_data = DD["corr_data"][:,:,:]

        if DD["flag"][0][0] == 1:
            flag = True
        else:
            flag = False

        if DD["py"][0][0] == 1:
            py = True
        else:
            py = False
        
        dt = DD["dt"][0][0]
        corr_method = DD["corr_method"][0]
        maxlag = DD["maxlag"][0][0]
        smoothspect_N = DD["smoothspect_N"][0][0]
        flag_gap = DD["flag_gap"][0][0]
        threads = DD["threads"][0][0]
        corr_pair = DD["corr_pair"][:,:]
        jobs = DD["jobs"][0][0]

        
    else:
        raise ValueError("format must be 'npz', 'h5' or 'mat'")
    
    CorrData = CorrData_Class(corr_data, dt, corr_method, corr_pair, maxlag, smoothspect_N,
        flag, flag_gap, threads, jobs, py)
    
    return CorrData



def load_stack(filename, format="npz", only_header=False):
    par = None

    if format == "npz":
        DD = np.load(filename, allow_pickle=True)
        if only_header:
            stack_data = None
        else:
            stack_data = DD["stack_data"]

        stack_ngood = DD["stack_ngood"]
        dt = DD["dt"].reshape(1)[0]
        stack_method = DD["stack_method"].reshape(1)[0]
        par = DD["par"].reshape(1)[0]
        stack_all = DD["stack_all"].reshape(1)[0]
        stack_len = DD["stack_len"].reshape(1)[0]
        stack_step = DD["stack_step"].reshape(1)[0]
        pick = DD["pick"].reshape(1)[0]
        median_high = DD["median_high"].reshape(1)[0]
        median_low = DD["median_low"].reshape(1)[0]
        flag = DD["flag"].reshape(1)[0]
        flag_gap = DD["flag_gap"].reshape(1)[0]
        threads = DD["threads"].reshape(1)[0]
        jobs = DD["jobs"].reshape(1)[0]
        py = DD["py"].reshape(1)[0]

    elif format == "h5":
        with h5py.File(filename, 'r') as f:
            group = f['noiseflow_group']

            if only_header:
                stack_data = None
            else:
                stack_data = group['stack_data'][:]

            dt = group.attrs['dt']
            stack_method = group.attrs['stack_method']
            stack_all = group.attrs['stack_all']
            stack_len = group.attrs['stack_len']
            stack_step = group.attrs['stack_step']
            pick = group.attrs['pick']
            median_high = group.attrs['median_high']
            median_low = group.attrs['median_low']
            flag = group.attrs['flag']
            flag_gap = group.attrs['flag_gap']
            threads = group.attrs['threads']
            stack_ngood = group['stack_ngood'][:]
            jobs = group.attrs['jobs']
            py = group.attrs['py']

            par = {}
            for key in group.attrs.keys():
                value = group.attrs[key]
                if isinstance(value, float) and np.isnan(value):
                    value = None
                par[key] = value

    elif format == "mat":
        DD = loadmat(filename)
        if only_header:
            stack_data = None
        else:
            stack_data = DD["stack_data"][:,:,:]

        if DD["stack_all"][0][0] == 1:
            stack_all = True
        else:
            stack_all = False

        if DD["pick"][0][0] == 1:
            pick = True
        else:
            pick = False

        if DD["flag"][0][0] == 1:
            flag = True
        else:
            flag = False
        
        if DD["py"][0][0] == 1:
            py = True
        else:
            py = False
        
        dt = DD["dt"][0][0]
        stack_method = DD["stack_method"][0]
        stack_len = DD["stack_len"][0][0]
        stack_step = DD["stack_step"][0][0]
        median_high = DD["median_high"][0][0]
        median_low = DD["median_low"][0][0]
        flag_gap = DD["flag_gap"][0][0]
        threads = DD["threads"][0][0]
        stack_ngood = DD["stack_ngood"][:]
        jobs = DD["jobs"][0][0]
        par = decode_dict_from_mat(DD['par'][0][0])

    else:
        raise ValueError("format must be 'npz', 'h5' or 'mat'")
    
    StackData = StackData_Class(stack_data, stack_ngood, dt,
                stack_method, par, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, flag_gap, threads, jobs, py)
    
    return StackData



def decode_dict_from_mat(par_data):
    par = {}

    for field_name in par_data.dtype.names:
        value = par_data[field_name][0, 0]
    
        if isinstance(value, np.ndarray) and value.size == 1:
            value = value.item()
    
        if isinstance(value, float) and np.isnan(value):
            value = None

        par[field_name] = value

    if par["stat"] == 1:
        par["stat"] = True
    else:
        par["stat"] = False

    if par["plot"] == 1:
        par["plot"] = True
    else:
        par["plot"] = False

    if par["normalize"] == 1:
        par["normalize"] = True
    else:
        par["normalize"] = False

    return par