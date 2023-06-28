import time
import scipy
import numpy as np

from tqdm import tqdm
from joblib import Parallel, delayed
from noiseflow.cc.python.utils import slice_window, moving_ave, whiten, split


def time_norm_clip(data, clip_std):
    lim = clip_std * np.std(data)
    time_norm_data = np.clip(data, -lim, lim)

    return time_norm_data


def time_norm_smooth(data, smooth_N):
    temp = moving_ave(np.abs(data), smooth_N)
    time_norm_data = np.divide(data, 
                               temp, 
                               out=np.zeros_like(data, dtype=data.dtype), 
                               where=temp!=0)

    return time_norm_data


class RFFTClass_python(object):
    def __init__(self, raw_data, dt, 
                cc_len, cc_step, 
                time_norm, clip_std, smooth_N,
                freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
                flag, jobs):
        
        self.raw_data = raw_data
        self.dt = dt
        self.cc_len = cc_len
        self.cc_step = cc_step
        self.time_norm_method = time_norm
        self.clip_std = clip_std
        self.smooth_N = smooth_N
        self.freq_norm_method = freq_norm
        self.freqmin = freqmin
        self.freqmax = freqmax
        self.whiten_npad = whiten_npad
        self.smoothspect_N = smoothspect_N
        self.flag = flag
        self.jobs = jobs

        # convert time to points
        self.cc_len_points = int(cc_len/dt)
        self.cc_step_points = int(cc_step/dt)
        self.npts = int(raw_data.shape[1])
        self.rfft_npts = int(self.cc_len_points/2+1)
        self.channel_num = int(raw_data.shape[0])
  
        # slice_window --> _win_num, _win_info
        self.win_info = slice_window(self.npts, self.cc_len_points, self.cc_step_points)
        self.win_num = int(self.win_info.shape[0])

        # initialize output_data
        if self.raw_data.dtype == np.dtype(np.float32):
            self.output_data = np.empty((self.channel_num, self.win_num, self.rfft_npts), dtype=np.complex64)
            self.time_dtype = np.float32
            self.freq_dtype = np.complex64
        elif self.raw_data.dtype == np.dtype(np.float64):
            self.output_data = np.empty((self.channel_num, self.win_num, self.rfft_npts), dtype=np.complex128)
            self.time_dtype = np.float64
            self.freq_dtype = np.complex128
        else:
            print("error: please input the correct data type.")
            exit(1)


    def time_norm(self, data):
        time_norm_data = None
        if self.time_norm_method == "no":
            time_norm_data = data
        elif self.time_norm_method == "clip":
            time_norm_data = time_norm_clip(data, self.clip_std)
        elif self.time_norm_method == "onebit":
            time_norm_data = np.sign(data)
        elif self.time_norm_method == "smooth":
            time_norm_data = time_norm_smooth(data, self.smooth_N)
        else:
            print("error: please input the correct time_norm method.")
            exit(1)

        return time_norm_data
    

    def rfft(self, time_norm_data):
        rfft_data = np.empty((self.win_num, self.rfft_npts), dtype=self.freq_dtype)

        for i in range(self.win_num):
            segment = time_norm_data[self.win_info[i, 0]:self.win_info[i, 1]]
            segment_rfft = scipy.fft.rfft(segment)
            rfft_data[i] = segment_rfft

        return rfft_data


    def freq_norm(self, rfft_data):
        freq_norm_data = None
        if self.freq_norm_method == "no":
            freq_norm_data = rfft_data
        else:
            freq_norm_data = whiten(rfft_data, self.dt, self.freq_norm_method, self.freqmin, self.freqmax, self.smoothspect_N, self.whiten_npad)

        return freq_norm_data


    def process_chunk(self, chunk_start, chunk_end):
        results = []
        if self.flag and chunk_start == 0:
            bar = tqdm(range(chunk_start, chunk_end))
        else:
            bar = range(chunk_start, chunk_end)

        for i in bar:
            # time_norm
            time_norm_data = self.time_norm(self.raw_data[i, :])

            # rfft
            rfft_data = self.rfft(time_norm_data)

            # freq_norm
            rfft_norm_data = self.freq_norm(rfft_data)

            results.append(rfft_norm_data)

        return results


    def run(self):
        if self.flag:
            start_time = time.time()
            print(f"Start rfft with {self.jobs} jobs in python...")
        
        # parallel processing
        if self.jobs > 1:
            chunk_list = split(num=self.channel_num, n_jobs=self.jobs)
            results = Parallel(n_jobs=self.jobs, backend="loky")(delayed(self.process_chunk)(chunk_start, chunk_end) for chunk_start, chunk_end in chunk_list)

            for i, chunk_results in enumerate(results):
                chunk_start = chunk_list[i][0]
                for j, result in enumerate(chunk_results):
                    self.output_data[chunk_start + j] = result
        else:
            results = self.process_chunk(0, self.channel_num)
            self.output_data = np.array(results, dtype=self.freq_dtype)

        if self.flag:
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"End rfft with total time {elapsed_time}s")

