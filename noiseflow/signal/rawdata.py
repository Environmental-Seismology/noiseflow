#%%
import h5py
import numpy as np
import matplotlib.pyplot as plt

from scipy.io import savemat
from noiseflow.signal.wrapper import bandpass, bandstop, lowpass, highpass, detrend, decimate, taper


class RawData_Class(object):
    def __init__(self, data, sampling_rate):
        self.data = data
        self.sampling_rate = sampling_rate


    def detrend(self, type='linear', flag=False, flag_gap=None, threads=1, py=False):
        detrend(self.data, type, flag, flag_gap, threads, py)


    def taper(self, max_percentage=0.05, type='hann', side='both', flag=False, flag_gap=None, threads=1, py=False):  
        taper(self.data, max_percentage, type, side, flag, flag_gap, threads, py)


    def bandpass(self, freqmin, freqmax, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
        bandpass(self.data, freqmin, freqmax, self.sampling_rate, corners, zerophase, flag, flag_gap, threads, py)


    def bandstop(self, freqmin, freqmax, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):  
        bandstop(self.data, freqmin, freqmax, self.sampling_rate, corners, zerophase, flag, flag_gap, threads, py)
        

    def lowpass(self, freq, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
        lowpass(self.data, freq, self.sampling_rate, corners, zerophase, flag, flag_gap, threads, py)
        

    def highpass(self, freq, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False):
        highpass(self.data, freq, self.sampling_rate, corners, zerophase, flag, flag_gap, threads, py)
        

    def decimate(self, sampling_rate_new, flag=False, flag_gap=None, threads=1, py=False):
        results = decimate(self.data, self.sampling_rate, sampling_rate_new, flag, flag_gap, threads, py)
        self.data = results
        self.sampling_rate = sampling_rate_new


    # save data
    def save(self, save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3):
        if format == 'npz':
            if compression:
                np.savez_compressed(save_path,
                                    data = self.data,
                                    sampling_rate = self.sampling_rate)
            else:
                np.savez(save_path,
                        data = self.data,
                        sampling_rate = self.sampling_rate)
                
        elif format == 'h5':
            with h5py.File(save_path, 'w') as f:
                group = f.create_group('noiseflow_group')
                group.attrs['sampling_rate'] = self.sampling_rate
                if compression:
                    group.create_dataset('data', data=self.data, compression=h5_compression_format, compression_opts=h5_compression_opts)
                else:
                    group.create_dataset('data', data=self.data)

        elif format == 'mat':
            savemat(save_path, 
                {'data': self.data,
                'sampling_rate': self.sampling_rate},
                do_compression=compression)
            
        else:
            raise ValueError('format must be npz, h5, or mat')


    def plot(self, channel_num=0, win_num=0, save=False, save_path=None, dpi=100):
        pass


