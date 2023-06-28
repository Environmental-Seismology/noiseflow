import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

from scipy.io import savemat
from obspy.imaging.spectrogram import spectrogram


class RFFTData_Class(object):
    def __init__(self,rfft_data, dt, cc_len, cc_step, 
                time_norm, clip_std, smooth_N,
                freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N,
                flag, flag_gap, threads, jobs, py):
        self.rfft_data = rfft_data
        self.dt = dt
        self.cc_len = cc_len
        self.cc_step=cc_step 
        self.time_norm=time_norm
        self.clip_std=clip_std
        self.smooth_N=smooth_N
        self.freq_norm=freq_norm
        self.freqmin=freqmin
        self.freqmax=freqmax
        self.whiten_npad=whiten_npad
        self.smoothspect_N=smoothspect_N
        self.flag=flag
        self.flag_gap=flag_gap
        self.threads=threads
        self.jobs=jobs
        self.py=py

    # save data
    def save(self, save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3):
        if format == 'npz':
            if compression:
                np.savez_compressed(save_path,
                        rfft_data=self.rfft_data,
                        dt=self.dt,
                        cc_len=self.cc_len,
                        cc_step=self.cc_step,
                        time_norm=self.time_norm,
                        clip_std=self.clip_std,
                        smooth_N=self.smooth_N,
                        freq_norm=self.freq_norm,
                        freqmin=self.freqmin,
                        freqmax=self.freqmax,
                        whiten_npad=self.whiten_npad,
                        smoothspect_N=self.smoothspect_N,
                        flag=self.flag,
                        flag_gap=self.flag_gap,
                        threads=self.threads,
                        jobs=self.jobs,
                        py=self.py)
            else:
                np.savez(save_path,
                        rfft_data=self.rfft_data,
                        dt=self.dt,
                        cc_len=self.cc_len,
                        cc_step=self.cc_step,
                        time_norm=self.time_norm,
                        clip_std=self.clip_std,
                        smooth_N=self.smooth_N,
                        freq_norm=self.freq_norm,
                        freqmin=self.freqmin,
                        freqmax=self.freqmax,
                        whiten_npad=self.whiten_npad,
                        smoothspect_N=self.smoothspect_N,
                        flag=self.flag,
                        flag_gap=self.flag_gap,
                        threads=self.threads,
                        jobs=self.jobs,
                        py=self.py)
                
        elif format == 'h5':
            with h5py.File(save_path, 'w') as f:
                group = f.create_group('noiseflow_group')
                group.attrs['dt'] = self.dt
                group.attrs['cc_len'] = self.cc_len
                group.attrs['cc_step'] = self.cc_step
                group.attrs['time_norm'] = self.time_norm
                group.attrs['clip_std'] = self.clip_std
                group.attrs['smooth_N'] = self.smooth_N
                group.attrs['freq_norm'] = self.freq_norm
                group.attrs['freqmin'] = self.freqmin
                group.attrs['freqmax'] = self.freqmax
                group.attrs['whiten_npad'] = self.whiten_npad
                group.attrs['smoothspect_N'] = self.smoothspect_N
                group.attrs['flag'] = self.flag
                group.attrs['flag_gap'] = self.flag_gap
                group.attrs['threads'] = self.threads
                group.attrs['jobs'] = self.jobs
                group.attrs['py'] = self.py
                if compression:
                    group.create_dataset('rfft_data', data=self.rfft_data, compression=h5_compression_format, compression_opts=h5_compression_opts)
                else:
                    group.create_dataset('rfft_data', data=self.rfft_data)

        elif format == 'mat':
            savemat(save_path, 
                {'rfft_data': self.rfft_data,
                'dt': self.dt,
                'cc_len': self.cc_len,
                'cc_step': self.cc_step,
                'time_norm': self.time_norm,
                'clip_std': self.clip_std,
                'smooth_N': self.smooth_N,
                'freq_norm': self.freq_norm,
                'freqmin': self.freqmin,
                'freqmax': self.freqmax,
                'whiten_npad': self.whiten_npad,
                'smoothspect_N': self.smoothspect_N,
                'flag': self.flag,
                'flag_gap': self.flag_gap,
                'threads': self.threads,
                'jobs': self.jobs,
                'py': self.py},
                do_compression=compression)
            
        else:
            raise ValueError('format must be npz, h5, or mat')


    # plot spectrum 
    def spectrogram(self, 
                    channel_indx=0, 
                    win_indx=0, 
                    raw_data=None,
                    dbscale=False,
                    log=True, 
                    figsize=(10, 4),
                    save=False, 
                    save_path=None, 
                    dpi=100):
        
        if raw_data is None:
            win1 = win_indx * (int(self.cc_len/self.dt) - int(self.cc_step/self.dt))
            win2 = win1 + int(self.cc_len/self.dt)
            rfft_whitedata = np.fft.irfft(self.rfft_data[channel_indx,win_indx]).real

            fig, ax = plt.subplots(figsize=figsize)
            spectrogram(rfft_whitedata, 1/self.dt,  axes=ax, dbscale=dbscale, log=log)

            ax.set_title("RFFTData: channel_indx=%d, win_indx=%d, time_norm=%s \nfreq_norm=%s, freq_band=[%.2f, %.2f] hz" % (channel_indx, win_indx, self.time_norm, self.freq_norm, self.freqmin, self.freqmax))
            ax.set_xlim(0, self.cc_len)
            ax.set_ylabel("Time(s)")
            ax.set_ylim(self.freqmin/2, self.freqmax*2)
            ax.set_ylabel("Frequency(hz)")
            ax.plot([0, self.cc_len], [self.freqmin, self.freqmin], '--', color='red', lw=3, alpha=0.7)
            ax.plot([0, self.cc_len], [self.freqmax, self.freqmax], '--', color='red', lw=3, alpha=0.7)

        else:
            win1 = win_indx * (int(self.cc_len/self.dt) - int(self.cc_step/self.dt))
            win2 = win1 + int(self.cc_len/self.dt)
            rfft_whitedata = np.fft.irfft(self.rfft_data[channel_indx,win_indx]).real

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
            spectrogram(raw_data[win1:win2], 1/self.dt, axes=ax1, dbscale=dbscale, log=log)
            spectrogram(rfft_whitedata, 1/self.dt, axes=ax2, dbscale=dbscale, log=log)

            fig.suptitle("RFFTData: channel_indx=%d, win_indx=%d, time_norm=%s \nfreq_norm=%s, freq_band=[%.2f, %.2f] hz" % (channel_indx, win_indx, self.time_norm, self.freq_norm, self.freqmin, self.freqmax))
            ax1.set_ylim(self.freqmin/2, self.freqmax*2)
            ax1.set_xlim(0, self.cc_len)
            ax1.set_ylabel("Frequency(hz)")
            ax1.plot([0, self.cc_len], [self.freqmin, self.freqmin], '--', color='red', lw=3, alpha=0.7)
            ax1.plot([0, self.cc_len], [self.freqmax, self.freqmax], '--', color='red', lw=3, alpha=0.7)
            ax1.set_title("no whitening")

            ax2.set_ylim(self.freqmin/2, self.freqmax*2)
            ax2.set_xlim(0, self.cc_len)
            ax2.set_ylabel("Frequency(hz)")
            ax2.set_xlabel("Time(s)")
            ax2.plot([0, self.cc_len], [self.freqmin, self.freqmin], '--', color='red', lw=3, alpha=0.7)
            ax2.plot([0, self.cc_len], [self.freqmax, self.freqmax], '--', color='red', lw=3, alpha=0.7)
            ax2.set_title("whitening")
        
        if save:
            if save_path is None:
                raise ValueError("save_path must be specified")
            fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()


