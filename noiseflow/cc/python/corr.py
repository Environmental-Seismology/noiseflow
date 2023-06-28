import time
import scipy
import numpy as np

from tqdm import tqdm 
from joblib import Parallel, delayed
from noiseflow.cc.python.utils import moving_ave, split


def xcorr(source_data, receiver_data):
    rfft_corr_data = np.conj(source_data) * receiver_data

    return rfft_corr_data


def deconv(source_data, receiver_data, smoothspect_N):
    win_num = int(source_data.shape[0])
    rfft_npts = int(source_data.shape[1])
    ss_data = np.empty((win_num, rfft_npts), dtype=source_data.dtype)

    for i in range(win_num):
        ss = source_data[i, :]
        temp = moving_ave(np.abs(ss), smoothspect_N)
        ss_data[i] = np.divide(ss, 
                               temp, 
                               out=np.zeros_like(ss, dtype=ss.dtype), 
                               where=temp!=0)

    rfft_corr_data = np.conj(ss_data) * receiver_data

    return rfft_corr_data


def coherency(source_data, receiver_data, smoothspect_N):
    win_num = int(source_data.shape[0])
    rfft_npts = int(source_data.shape[1])

    ss_data = np.empty((win_num, rfft_npts), dtype=source_data.dtype)
    rr_data = np.empty((win_num, rfft_npts), dtype=source_data.dtype)

    for i in range(win_num):
        ss = source_data[i, :]
        rr = receiver_data[i, :]
        ss_temp = moving_ave(np.abs(ss), smoothspect_N)
        rr_temp = moving_ave(np.abs(rr), smoothspect_N)
        ss_data[i] = np.divide(ss, 
                                ss_temp, 
                                out=np.zeros_like(ss, dtype=ss.dtype), 
                                where=ss_temp!=0)
        
        rr_data[i] = np.divide(rr,
                                rr_temp,
                                out=np.zeros_like(rr, dtype=rr.dtype),
                                where=rr_temp!=0)


    rfft_corr_data = np.conj(ss_data) * rr_data

    return rfft_corr_data


class CorrClass_python(object):
    def __init__(self, rfft_data, dt, corr_method, corr_pair, maxlag,
                smoothspect_N, flag, jobs):
        
        self.rfft_data = rfft_data
        self.dt = dt
        self.corr_method = corr_method
        self.corr_pair = corr_pair
        self.maxlag = maxlag
        self.smoothspect_N = smoothspect_N
        self.flag = flag
        self.jobs = jobs

        self.pair_num = int(corr_pair.shape[0])
        self.channel_num = int(rfft_data.shape[0])
        self.win_num = int(rfft_data.shape[1])
        self.rfft_npts = int(rfft_data.shape[2])
        self.npts = int(2*(self.rfft_npts-1))

        t = np.arange(-self.rfft_npts + 1, self.rfft_npts) * dt
        ind = np.where(np.abs(t) <= maxlag)[0]
        self.maxlag_npts = int(ind.size)
        self.maxlag_start = int(ind[0])
        self.maxlag_end = int(ind[-1] + 1)

        # initialize output_data
        if self.rfft_data.dtype == np.dtype(np.complex64):
            self.output_data = np.empty((self.pair_num, self.win_num, self.maxlag_npts), dtype=np.float32)
            self.time_dtype = np.float32
            self.freq_dtype = np.complex64
        elif self.rfft_data.dtype == np.dtype(np.complex128):
            self.output_data = np.empty((self.pair_num, self.win_num, self.maxlag_npts), dtype=np.float64)
            self.time_dtype = np.float64
            self.freq_dtype = np.complex128
        else:
            print("error: please input the correct data type.")
            exit(1)


    def corr(self, source_data, receiver_data):
        if self.corr_method == "xcorr":
            rfft_corr_data = xcorr(source_data, receiver_data)
        elif self.corr_method == "deconv":
            rfft_corr_data = deconv(source_data, receiver_data, self.smoothspect_N)
        elif self.corr_method == "coherency":
            rfft_corr_data = coherency(source_data, receiver_data, self.smoothspect_N)
        else:
            print("error: please input the correct corr method.")
            exit(1)
        
        return rfft_corr_data


    def irfft(self, rfft_corr_data):
        corr_data = np.empty((self.win_num, self.npts), dtype=self.time_dtype)


        for i in range(self.win_num):
            ss = rfft_corr_data[i]
            corr_data[i] = np.roll(scipy.fft.irfft(ss), int(self.rfft_npts - 1))

        return corr_data


    def cut_maxlag(self, corr_data):
        maxlag_corr_data = corr_data[:, self.maxlag_start:self.maxlag_end]

        return maxlag_corr_data


    def process_chunk(self, chunk_start, chunk_end):
        results = []
        if self.flag and chunk_start == 0:
            bar = tqdm(range(chunk_start, chunk_end))
        else:
            bar = range(chunk_start, chunk_end)

        for i in bar:
            # corr
            source_data = self.rfft_data[self.corr_pair[i, 0], :, :]
            receiver_data = self.rfft_data[self.corr_pair[i, 1], :, :]
            rfft_corr_data = self.corr(source_data, receiver_data)

            # irfft
            corr_data = self.irfft(rfft_corr_data)

            # cut maxlag
            maxlag_corr_data = self.cut_maxlag(corr_data)

            results.append(maxlag_corr_data)

        return results


    def run(self):
        if self.flag:
            start_time = time.time()
            print(f"Start corr with {self.jobs} jobs in python...")
        
        # parallel processing
        if self.jobs > 1:
            chunk_list = split(num=self.pair_num, n_jobs=self.jobs)
            results = Parallel(n_jobs=self.jobs, backend="loky")(delayed(self.process_chunk)(chunk_start, chunk_end) for chunk_start, chunk_end in chunk_list)

            for i, chunk_results in enumerate(results):
                chunk_start = chunk_list[i][0]
                for j, result in enumerate(chunk_results):
                    self.output_data[chunk_start + j] = result
        else:
            results = self.process_chunk(0, self.pair_num)
            self.output_data = np.array(results, dtype=self.time_dtype)

        if self.flag:
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"End corr with total time {elapsed_time}s")



