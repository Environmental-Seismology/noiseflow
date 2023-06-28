import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

from scipy.io import savemat
from obspy.core import UTCDateTime
from obspy.signal.filter import bandpass
from noiseflow.cc.utils_time import time_linspace, get_timestamp


class CorrData_Class(object):
    def __init__(self, corr_data, dt, 
                 corr_method, corr_pair, maxlag,smoothspect_N,
                flag, flag_gap, threads, jobs, py):
        self.corr_data = corr_data
        self.dt = dt
        self.corr_method = corr_method
        self.corr_pair = corr_pair
        self.maxlag = maxlag
        self.smoothspect_N=smoothspect_N
        self.flag=flag
        self.flag_gap=flag_gap
        self.threads=threads
        self.jobs=jobs
        self.py=py

    def save(self, save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3):
        if format == 'npz':
            if compression:
                np.savez_compressed(save_path,
                                    corr_data=self.corr_data,
                                    dt=self.dt,
                                    corr_method = self.corr_method,
                                    corr_pair = self.corr_pair,
                                    maxlag = self.maxlag,
                                    smoothspect_N = self.smoothspect_N,
                                    flag = self.flag,
                                    flag_gap = self.flag_gap,
                                    threads = self.threads,
                                    jobs = self.jobs,
                                    py = self.py)
            else:
                np.savez(save_path,
                        corr_data=self.corr_data,
                        dt=self.dt,
                        corr_method = self.corr_method,
                        corr_pair = self.corr_pair,
                        maxlag = self.maxlag,
                        smoothspect_N = self.smoothspect_N,
                        flag = self.flag,
                        flag_gap = self.flag_gap,
                        threads = self.threads,
                        jobs = self.jobs,
                        py = self.py)
                
        elif format == 'h5':
            with h5py.File(save_path, 'w') as f:
                group = f.create_group('noiseflow_group')
                group.attrs['dt'] = self.dt
                group.attrs['corr_method'] = self.corr_method
                group.attrs['corr_pair'] = self.corr_pair
                group.attrs['maxlag'] = self.maxlag
                group.attrs['smoothspect_N'] = self.smoothspect_N
                group.attrs['flag'] = self.flag
                group.attrs['flag_gap'] = self.flag_gap
                group.attrs['threads'] = self.threads
                group.attrs['jobs'] = self.jobs
                group.attrs['py'] = self.py
                group.create_dataset('corr_pair', data=self.corr_pair)
                if compression:
                    group.create_dataset('corr_data', data=self.corr_data, compression=h5_compression_format, compression_opts=h5_compression_opts)
                else:
                    group.create_dataset('corr_data', data=self.corr_data)

        elif format == 'mat':
            savemat(save_path,
                {'corr_data': self.corr_data,
                 'dt': self.dt,
                 'corr_method': self.corr_method,
                 'corr_pair': self.corr_pair,
                 'maxlag': self.maxlag,
                 'smoothspect_N': self.smoothspect_N,
                 'flag': self.flag,
                 'flag_gap': self.flag_gap,
                 'threads': self.threads,
                 'jobs': self.jobs,
                 'py': self.py},
                do_compression=compression) 
        else:
            raise ValueError("format must be 'npz', 'h5' or 'mat'")
        

    # plot
    def plot(self, 
             pair_indx=0, 
             t_min=UTCDateTime("1970-01-01T00:00:00.0"), 
             cc_len=None, 
             cc_step=None, 
             win_start=None, 
             win_end=None, 
             lag_start=None, 
             lag_end=None, 
             amp_normalize=True, 
             amp_scale=1, 
             filter=False, 
             f1=None, 
             f2=None, 
             corners=4, 
             zerophase=True, 
             win_interval=None, 
             mode='waveform', 
             cmap='seismic',
             linewidth=0.8, 
             yticklabel_num=5, 
             figsize=(10, 6),
             save=False, 
             save_path=None, 
             dpi=100):
        
        # check pair_indx
        if pair_indx >= self.corr_data.shape[0]:
            raise ValueError("pair_indx must be <= (corr_data.shape[0]=%d)" % self.corr_data.shape[0])
        
        # init cc_len cc_step
        if cc_len == None:
            cc_len = self.corr_data.shape[2]*self.dt
        if cc_step == None:
            cc_step = 0.0

        # init timestamp
        t_max, win_interval_raw, win_vector = get_timestamp(win_num=self.corr_data.shape[1], 
                                                        cc_len=cc_len, 
                                                        cc_step=cc_step, 
                                                        t_min=t_min)

        # init win_start
        if win_start == None:
            win_start = t_min
        if win_start < t_min:
            raise ValueError("win_start must be >= (t_min=%s)" % str(t_min))

        # init win_end
        if win_end == None:
            win_end = t_max
        if win_end > t_max:
            raise ValueError("win_end must be <= (t_max=%s)" % str(t_max))

        # init win_interval
        if win_interval == None:
            win_interval = win_interval_raw      
        if win_interval < win_interval_raw:
            raise ValueError("win_interval must be >= %f" % win_interval_raw)
        if win_interval > t_max-t_min-cc_len:
            raise ValueError("win_interval must be <= %f" % (t_max-t_min-cc_len))

        # init lag_start 
        if lag_start == None:
            lag_start = -self.maxlag
        if lag_start < -self.maxlag:
            raise ValueError("lag_start must be >= (-maxlag=-%f)" % self.maxlag)
        
        # init lag_end
        if lag_end == None:
            lag_end = self.maxlag
        if lag_end > self.maxlag:
            raise ValueError("lag_end must be <= (maxlag=%f)" % self.maxlag)
        
        # lag vector
        tt = np.arange(lag_start, lag_end+self.dt, self.dt)
        tp_start = round((lag_start+self.maxlag)/self.dt)
        tp_end = tp_start+len(tt)  # make sure same length, so do not use round((lag_end+self.maxlag+self.dt)/self.dt)

        # select and filter data according to win_start, win_end
        win_indx = np.where((win_vector >= win_start) & (win_vector <= win_end))[0]
        if filter:
            data=np.empty((len(win_indx), len(tt)))
            for i in range(0, len(win_indx)):
                data[i,:] = bandpass(self.corr_data[pair_indx, win_indx[i], tp_start:tp_end], f1, f2, int(1/self.dt), corners=corners, zerophase=zerophase)
        else:
            data = self.corr_data[pair_indx, win_indx, tp_start:tp_end]

        # plot
        if mode == 'waveform':
            # init fig
            fig, ax = plt.subplots(figsize=figsize)

            # normalize
            if amp_normalize:
                scale = np.max(np.abs(data), axis=1)[:,None]
                if np.min(scale) != 0:
                    data_plot = amp_scale*data/scale
                else:
                    raise ValueError("data is all zeros in amp_normalize=True")
            else:
                scale = np.max(np.abs(data))
                if scale != 0:
                    data_plot = amp_scale*data/scale
                else:
                    raise ValueError("data is all zeros in amp_normalize=False")

            # get win_interval index
            ydist_n=[]
            win_interval_n = round(win_interval/win_interval_raw)
            for i in range(0, len(win_indx), win_interval_n):
                ydist_n.append(i)
            ydist_n = np.int64(np.array(ydist_n))

            # plot wavefrom
            for i in range(0, len(ydist_n)):
                ax.plot(tt, data_plot[ydist_n[i]]+ydist_n[i], 'k', linewidth=linewidth)

            # set x axis
            ax.set_title("CorrData: pair_indx=%d, filter=[%.2f, %.2f] hz" % (pair_indx, f1, f2))
            ax.set_xlim(tt[0], tt[-1])
            ax.set_xlabel('Time(s)')
            
            # set y axis
            yy_tick = np.linspace(0, len(win_indx)-1, num=yticklabel_num)
            yy_label_UTC = time_linspace(win_start, win_end, num=yticklabel_num)
            yy_label = np.array([i.strftime('%Y-%m-%dT%H:%M:%S') for i in yy_label_UTC])
            ax.set_yticks(yy_tick)
            ax.set_yticklabels(yy_label, rotation=45)

            # save or show
            if save:
                if save_path is None:
                    raise ValueError("save_path must be specified")
                fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches='tight')
                plt.close(fig)
            else:
                plt.show()

        elif mode == 'mat':
            # init fig
            fig, ax = plt.subplots(figsize=figsize)

            # average all data
            cc_mean = np.mean(data, axis=0)/np.max(np.mean(data, axis=0))

            # average win_interval data
            ntrace = int((win_end-win_start)/win_interval)
            ndata  = np.zeros(shape=(ntrace,len(tt)))
            for i in range(0,ntrace):
                tindx = np.where(((win_vector-win_start)>=i*win_interval) & ((win_vector-win_start)<(i+1)*win_interval))[0]
                ndata[i] = np.mean(data[tindx],axis=0)

            # normalize waveforms
            if amp_normalize:
                scale = np.max(np.abs(ndata), axis=1)[:,None]
                if np.min(scale) != 0:
                    ndata_plot = ndata/scale
                else:
                    raise ValueError("data is all zeros in amp_normalize=True")
            else:
                scale = np.max(np.abs(ndata))
                if scale != 0:
                    ndata_plot = ndata/scale
                else:
                    raise ValueError("data is all zeros in amp_normalize=False")
                
            # plot mat
            cax=ax.matshow(ndata_plot, extent=[lag_start, lag_end, ntrace, 0], aspect='auto', cmap=cmap)
            ax.plot(tt, ntrace/10*cc_mean+ntrace/2, 'k', linewidth=linewidth)
            
            # set x axis
            ax.set_title("CorrData: pair_indx=%d, filter=[%.2f, %.2f] hz" % (pair_indx, f1, f2))
            ax.set_xlim(tt[0], tt[-1])
            ax.set_xlabel('Time(s)')
            ax.xaxis.set_ticks_position('bottom')
            
            # set y axis
            ax.invert_yaxis()
            ax.set_ylim(0, ntrace)
            yy_tick=np.linspace(0, ntrace, num=yticklabel_num)
            yy_label_UTC=time_linspace(win_start, win_end, num=yticklabel_num)
            yy_label = np.array([i.strftime('%Y-%m-%dT%H:%M:%S') for i in yy_label_UTC])
            ax.set_yticks(yy_tick)
            ax.set_yticklabels(yy_label, rotation=45)

            # save or show
            if save:
                if save_path is None:
                    raise ValueError("save_path must be specified")
                fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches = 'tight') 
                plt.close(fig)
            else:
                plt.show()

        else:
            raise ValueError("mode must be 'waveform' or 'mat'")
        
