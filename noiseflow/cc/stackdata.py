import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from scipy.io import savemat
from obspy.core import UTCDateTime
from obspy.signal.filter import bandpass
from noiseflow.cc.utils_time import time_linspace


class StackData_Class(object):
    def __init__(self, stack_data, stack_ngood, dt,
                 stack_method, par, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, flag_gap, threads, jobs, py):
        self.stack_data = stack_data
        self.stack_ngood = stack_ngood
        self.dt = dt
        self.stack_method = stack_method
        self.par = par
        self.stack_all = stack_all
        self.stack_len = stack_len
        self.stack_step = stack_step
        self.pick = pick
        self.median_high = median_high
        self.median_low = median_low
        self.flag=flag
        self.flag_gap=flag_gap
        self.threads=threads
        self.jobs=jobs
        self.py=py

    # save data
    def save(self, save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3):
        self.maxlag = self.dt * (self.stack_data.shape[2]-1)/2

        if format == 'npz':
            if compression:
                np.savez_compressed(save_path,
                                    stack_data = self.stack_data,
                                    stack_ngood = self.stack_ngood,
                                    dt = self.dt,
                                    stack_method = self.stack_method,
                                    par = self.par,
                                    stack_all = self.stack_all,
                                    stack_len = self.stack_len,
                                    stack_step = self.stack_step,
                                    pick = self.pick,
                                    median_high = self.median_high,
                                    median_low = self.median_low,
                                    flag = self.flag,
                                    flag_gap = self.flag_gap,
                                    threads = self.threads,
                                    jobs = self.jobs,
                                    py = self.py,
                                    maxlag = self.maxlag)
            else:
                np.savez(save_path,
                        stack_data = self.stack_data,
                        stack_ngood = self.stack_ngood,
                        dt = self.dt,
                        stack_method = self.stack_method,
                        par = self.par,
                        stack_all = self.stack_all,
                        stack_len = self.stack_len,
                        stack_step = self.stack_step,
                        pick = self.pick,
                        median_high = self.median_high,
                        median_low = self.median_low,
                        flag = self.flag,
                        flag_gap = self.flag_gap,
                        threads = self.threads,
                        jobs = self.jobs,
                        py = self.py,
                        maxlag = self.maxlag)
                
        elif format == 'h5':
            par_copy = self.par.copy()
            with h5py.File(save_path, 'w') as f:
                group = f.create_group('noiseflow_group')
                group.attrs['dt'] = self.dt
                group.attrs['stack_method'] = self.stack_method
                group.attrs['stack_all'] = self.stack_all
                group.attrs['stack_len'] = self.stack_len
                group.attrs['stack_step'] = self.stack_step
                group.attrs['pick'] = self.pick
                group.attrs['median_high'] = self.median_high
                group.attrs['median_low'] = self.median_low
                group.attrs['flag'] = self.flag
                group.attrs['flag_gap'] = self.flag_gap
                group.attrs['threads'] = self.threads
                group.attrs['jobs'] = self.jobs
                group.attrs['py'] = self.py
                group.attrs['maxlag'] = self.maxlag

                for key, value in par_copy.items():
                    if value is None:
                        value = float('nan')
                    group.attrs[key] = value

                group.create_dataset('stack_ngood', data=self.stack_ngood)
                if compression:
                    group.create_dataset('stack_data', data=self.stack_data, compression=h5_compression_format, compression_opts=h5_compression_opts)
                else:
                    group.create_dataset('stack_data', data=self.stack_data)

        elif format == 'mat':
            par_copy = self.par.copy()
            for key, value in par_copy.items():
                if value is None:
                    par_copy[key] = np.nan

            savemat(save_path,
                {'stack_data': self.stack_data,
                'stack_ngood': self.stack_ngood,
                'dt': self.dt,
                'stack_method': self.stack_method,
                'par': par_copy,
                'stack_all': self.stack_all,
                'stack_len': self.stack_len,
                'stack_step': self.stack_step,
                'pick': self.pick,
                'median_high': self.median_high,
                'median_low': self.median_low,
                'flag': self.flag,
                'flag_gap': self.flag_gap,
                'threads': self.threads,
                'jobs': self.jobs,
                'py': self.py,
                'maxlag': self.maxlag},
                do_compression=compression)
        else:
            raise ValueError("format must be 'npz', 'h5' or 'mat'")
        

    # plot
    def plot(self, 
             pair_indx=0, 
             t_min=UTCDateTime("1970-01-01T00:00:00.0"), 
             stack_len=1, 
             stack_step=0, 
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
             ngood_label=False, 
             save=False, 
             save_path=None, 
             dpi=300): 
        
        self.maxlag = self.dt * (self.stack_data.shape[2]-1)/2

        # check pair_indx
        if pair_indx >= self.stack_data.shape[0]:
            raise ValueError("pair_indx must be smaller than stack_data.shape[0]=%d" % self.stack_data.shape[0])

        # init cc_len cc_step
        if cc_len == None:
            cc_len = self.stack_data.shape[2]*self.dt
        if cc_step == None:
            cc_step = 0.0

        # define time vector
        stack_win_num = self.stack_data.shape[1]
        corr_win_num = (stack_len-stack_step)*(stack_win_num-1) + stack_len
        t_max = t_min + (corr_win_num-1)*(cc_len-cc_step) + cc_len
        stack_interval = (stack_len-stack_step)*(cc_len-cc_step)
        win_vector = np.array([t_min + i*stack_interval for i in range(0, stack_win_num)])

        # init win_start
        if win_start == None:
            win_start = t_min
        if win_start < t_min:
            raise ValueError("win_start must be larger than t_min=%s" % str(t_min))

        # init win_end
        if win_end == None:
            win_end = t_max
        if win_end > t_max:
            raise ValueError("win_end must be smaller than t_max=%s" % str(t_max))

        # init lag_start 
        if lag_start == None:
            lag_start = -self.maxlag
        if lag_start < -self.maxlag:
            raise ValueError("lag_start must be larger than -maxlag=-%f" % self.maxlag)
        
        # init lag_end
        if lag_end == None:
            lag_end = self.maxlag
        if lag_end > self.maxlag:
            raise ValueError("lag_end must be smaller than maxlag=%f" % self.maxlag)
        
        # lag vector
        tt = np.arange(lag_start, lag_end+self.dt, self.dt)
        tp_start = round((lag_start+self.maxlag)/self.dt)
        tp_end = tp_start+len(tt)  # make sure same length, so do not use round((lag_end+self.maxlag+self.dt)/self.dt)

        # check win_interval
        if win_interval == None:
            win_interval = (stack_len-stack_step)*(cc_len-cc_step)
        if win_interval < (stack_len-stack_step)*(cc_len-cc_step):
            raise ValueError("win_interval must be larger than (stack_len-stack_step)*(cc_len-cc_step)=%f" % (stack_len-stack_step)*(cc_len-cc_step))
        if win_interval > t_max-stack_len*(cc_len-cc_step)-t_min:
            raise ValueError("win_interval must be smaller than t_max-stack_len*(cc_len-cc_step)-t_min=%f" % (t_max-stack_len*(cc_len-cc_step)-t_min))
        
        # select win_indx data and filter --> data, ngood
        win_indx = np.where((win_vector >= win_start) & (win_vector <= win_end))[0]
        ngood = self.stack_ngood[pair_indx, win_indx]
        if filter:
            data=np.empty((len(win_indx), len(tt)))
            for i in range(0, len(win_indx)):
                data[i,:] = bandpass(self.stack_data[pair_indx, win_indx[i], tp_start:tp_end], f1, f2, int(1/self.dt), corners=corners, zerophase=zerophase)
        else:
            data = self.stack_data[pair_indx, win_indx, tp_start:tp_end]

        # plot
        if mode == 'waveform':
            # init figure
            fig, ax = plt.subplots(1, 2, gridspec_kw=dict(width_ratios=[5, 1]), sharey="row", figsize=figsize)
            fig.subplots_adjust(wspace=0.01)

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

            # select win_interval data
            ydist_n=[]
            ngood_n=[]
            win_interval_n = round(win_interval/((stack_len-stack_step)*(cc_len-cc_step)))
            for i in range(0, len(win_indx), win_interval_n):
                ydist_n.append(i)
                ngood_n.append(ngood[i])
            ydist_n = np.array(ydist_n)
            ngood_n = np.array(ngood_n)

            # plot wavefrom
            for i in range(0, len(ydist_n)):
                ax[0].plot(tt, data_plot[int(ydist_n[i])]+ydist_n[i], 'k', linewidth=linewidth)
            
            # set x axis
            ax[0].set_title("StackData: pair_indx=%d, filter=[%.2f, %.2f] hz" % (pair_indx, f1, f2))
            ax[0].set_xlim(tt[0], tt[-1])
            ax[0].set_xlabel('Time(s)')
            
            # set y axis
            yy_tick = np.linspace(0, len(win_indx)-1, num=yticklabel_num)
            yy_label_UTC = time_linspace(win_start, win_end, num=yticklabel_num)
            yy_label = np.array([i.strftime('%Y-%m-%dT%H:%M:%S') for i in yy_label_UTC])
            ax[0].set_yticks(yy_tick)
            ax[0].set_yticklabels(yy_label, rotation=45)

            # plot ngood 
            (markers, stemlines, baseline) = ax[1].stem(ydist_n, ngood_n, orientation='horizontal')
            plt.setp(markers, marker='D', markeredgecolor="orange", markeredgewidth=1.5)
            plt.setp(stemlines, linestyle="-", color="olive", linewidth=linewidth)
            if ngood_label:
                for i in range(0, len(ngood_n)):
                    ax[1].text(ngood_n[i]+np.max(ngood_n)/10, ydist_n[i], int(ngood_n[i]), va='center', ha='left')

            ax[1].set_xlim(0, np.max(ngood_n)*1.2)
            ax[1].set_xlabel('Ngood')
            ax[1].get_yaxis().set_visible(False)
            ax[1].spines['top'].set_visible(False)
            ax[1].spines['right'].set_visible(False)
            ax[1].spines['bottom'].set_visible(True)
            ax[1].spines['left'].set_visible(True)

            # save or show
            if save:
                if save_path is None:
                    raise ValueError("save_path must be specified")
                fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches='tight')
                plt.close(fig)
            else:
                plt.show()

        elif mode == 'mat':
            # init figure
            fig, ax = plt.subplots(1, 2, gridspec_kw=dict(width_ratios=[5, 1]), sharey="row", figsize=figsize)
            fig.subplots_adjust(wspace=0.01)

            # average waveforms
            cc_mean = np.mean(data, axis=0)/np.max(np.mean(data, axis=0))
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
                
            ### plot mat
            cax=ax[0].matshow(ndata_plot, extent=[lag_start, lag_end, ntrace, 0], aspect='auto', cmap=cmap)
            ax[0].plot(tt,ntrace/10*cc_mean+ntrace/2,'k',linewidth=linewidth)
            
            # set x axis
            ax[0].set_title("StackData: pair_indx=%d, filter=[%.2f, %.2f] hz" % (pair_indx, f1, f2))
            ax[0].set_xlim(tt[0], tt[-1])
            ax[0].set_xlabel('Time(s)')
            ax[0].xaxis.set_ticks_position('bottom')

            # set y axis
            ax[0].invert_yaxis()
            ax[0].set_ylim(0, ntrace)
            yy_tick = np.linspace(0, ntrace, num=yticklabel_num)
            yy_label_UTC = time_linspace(win_start, win_end, num=yticklabel_num)
            yy_label = np.array([i.strftime('%Y-%m-%dT%H:%M:%S') for i in yy_label_UTC])
            ax[0].set_yticks(yy_tick)
            ax[0].set_yticklabels(yy_label, rotation=45)

            ### plot ngood
            ngood = self.stack_ngood[pair_indx, win_indx]
            yy = np.linspace(0, ntrace, len(ngood))
            (markers, stemlines, baseline) = ax[1].stem(yy, ngood, orientation='horizontal')
            plt.setp(markers, marker='D', markeredgecolor="orange", markeredgewidth=1.5)
            plt.setp(stemlines, linestyle="-", color="olive", linewidth=linewidth)
            if ngood_label:
                for i in range(0, len(ngood)):
                    ax[1].text(ngood[i]+np.max(ngood)/10, yy[i], int(ngood[i]), va='center', ha='left')

            ax[1].set_xlim(0, np.max(ngood)*1.2)
            ax[1].set_xlabel('Ngood')
            ax[1].get_yaxis().set_visible(False)
            ax[1].spines['top'].set_visible(False)
            ax[1].spines['right'].set_visible(False)
            ax[1].spines['bottom'].set_visible(True)
            ax[1].spines['left'].set_visible(True)
            
            # save
            if save:
                if save_path is None:
                    raise ValueError("save_path must be specified")
                fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches = 'tight') 
                plt.close(fig)
            else:
                plt.show()

        else:
            raise ValueError("mode must be 'waveform' or 'mat'")
        


    def plot_moveout(self, 
                corr_pair, 
                pair_dist, 
                source_indx=None, 
                receiver_indx=None, 
                dist_start=None, 
                dist_end=None, 
                amp_scale=1, 
                amp_normalize=True,
                win_num = 0, 
                lag_start=None, 
                lag_end=None, 
                filter=False, 
                f1=None, 
                f2=None, 
                corners=4, 
                zerophase=True, 
                dist_interval=None, 
                mode='waveform', 
                cmap='seismic',
                linewidth=0.8, 
                yticklabel_num=10, 
                figsize=(10, 6),
                dist_unit="m", 
                velocity=[], 
                save=False, 
                save_path=None, 
                dpi=100): 
        
        self.maxlag = self.dt * (self.stack_data.shape[2]-1)/2

        # check source_indx and receiver_indx
        if source_indx != None and receiver_indx != None:
            raise ValueError("source_indx and receiver_indx cannot be specified at the same time")
        elif source_indx != None:
            if source_indx not in corr_pair[:,0]:
                raise ValueError("source_indx must be in corr_pair[:,0]")
        elif receiver_indx != None:
            if receiver_indx not in corr_pair[:,1]:
                raise ValueError("receiver_indx must be in corr_pair[:,1]")
        else:
            raise ValueError("source_indx or receiver_indx must be specified")
        
        # check corr_pair pair_dist
        if corr_pair.shape[0] != len(pair_dist):
            raise ValueError("the length of corr_pair is not consistent with pair_dist")
        
        # check demision of pair_dist
        if pair_dist.ndim == 2:
            pair_dist = pair_dist.reshape(-1)

        # check dist_unit
        if dist_unit not in ['m', 'km', 'degree']:
            raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
        
        # init dist_start
        if dist_start == None:
            dist_start = np.min(pair_dist)
        if dist_start < np.min(pair_dist):
            raise ValueError("dist_start must be >= (np.min(pair_dist)=%f)" % np.min(pair_dist))
        
        # init dist_end
        if dist_end == None:
            dist_end = np.max(pair_dist)
        if dist_end > np.max(pair_dist):
            raise ValueError("dist_end must be <= (np.max(pair_dist)=%f)" % np.max(pair_dist))

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

        # select dist_start dist_end
        if source_indx != None:
            a_id = np.where((corr_pair[:,0]==source_indx))[0]
            b_id = np.where((pair_dist>=dist_start) & (pair_dist<=dist_end))[0]
            id = np.intersect1d(a_id ,b_id)
        elif receiver_indx != None:
            a_id = np.where((corr_pair[:,1]==receiver_indx))[0]
            b_id = np.where((pair_dist>=dist_start) & (pair_dist<=dist_end))[0]
            id = np.intersect1d(a_id, b_id)
        else:
            raise ValueError("source_indx or receiver_indx must be specified")
        dist = pair_dist[id]

        # init dist_interval
        if dist_interval == None:
            dist_interval = np.min(dist)
        if dist_interval < np.diff(np.sort(dist)).min():
            raise ValueError("win_interval must be >= (np.diff(np.sort(dist)).min()=%f)" % np.diff(np.sort(dist)).min())
        if dist_interval > np.max(dist):
            raise ValueError("win_interval must be <= (np.max(dist)=%f)" % np.max(dist))
        
        # filter --> data
        if filter:
            data=np.empty((len(id), len(tt)))
            for i in range(0, len(id)):
                data[i,:] = bandpass(self.stack_data[id[i], win_num, tp_start:tp_end], f1, f2, int(1/self.dt), corners=corners, zerophase=zerophase)
        else:
            data = self.stack_data[id, win_num, tp_start:tp_end]
            
        # plot
        colors = list(mcolors.TABLEAU_COLORS.keys())
        if mode == 'waveform':
            if data.size != 0:
                # init fig
                fig, ax = plt.subplots(figsize=figsize)

                # normalize data
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
                
                # plot waveform
                for i in range(0, data_plot.shape[0]):
                    ax.plot(tt, amp_scale*data_plot[i]+dist[i], 'k', linewidth=linewidth)

                # plot velocity
                for i in range(0, len(velocity)):
                    x0=0; y0=0
                    if dist_unit == 'm':
                        x1 = (dist_end-dist_start)/velocity[i]
                        y1 = dist_end
                    elif dist_unit == 'km':
                        x1 = (dist_end-dist_start)/1000/velocity[i]
                        y1 = dist_end
                    elif dist_unit == 'degree':
                        x1 = (dist_end-dist_start)/(111.2*1000)/velocity[i]
                        y1 = dist_end
                    else:
                        raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
                    ax.plot([x0,x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5, label=str(velocity[i])+'m/s')
                    ax.plot([x0,-x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5)

                # set legend
                ax.legend(loc='upper right', fontsize=8, shadow=False)
                ax.set_xlim(tt[0], tt[-1])
                ax.set_xlabel('Time(s)')
                ax.set_ylabel('Distance(%s)' % (dist_unit))

                # set title
                if source_indx != None:
                    ax.set_title("StackData: source_indx=%d, win_num=%d, filter=[%.2f, %.2f] hz" % (source_indx, win_num, f1, f2))
                elif receiver_indx != None:
                    ax.set_title("StackData: receiver_indx=%d, win_num=%d, filter=[%.2f, %.2f] hz" % (source_indx, win_num, f1, f2))
                else:
                    raise ValueError("source_indx or receiver_indx must be specified")
            
                # save or show
                if save:
                    if save_path is None:
                        raise ValueError("save_path must be specified")
                    fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches = 'tight') 
                    plt.close(fig)
                else:
                    plt.show()
            else:
                raise ValueError("data is empty")

        elif mode=="mat":
            if data.size != 0:
                # init fig
                fig, ax = plt.subplots(figsize=figsize)

                # average waveforms
                ntrace = int((dist_end-dist_start)/dist_interval)
                ndata  = np.zeros(shape=(ntrace,len(tt)))
                for i in range(0,ntrace):
                    tindx = np.where(((dist-dist_start)>=i*dist_interval) & ((dist-dist_start)<(i+1)*dist_interval))[0]
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
                
                # plot
                cax=ax.matshow(ndata_plot, extent=[lag_start, lag_end, ntrace, 0], aspect='auto', cmap=cmap) # cmap=cmaps.MPL_jet,
                for i in range(0, len(velocity)):
                    x0=0
                    y0=0
                    if dist_unit == 'm':
                        x1 = (dist_end-dist_start)/velocity[i]
                        y1 = ntrace
                    elif dist_unit == 'km':
                        x1 = (dist_end-dist_start)/1000/velocity[i]
                        y1 = ntrace
                    elif dist_unit == 'degree':
                        x1 = (dist_end-dist_start)/(111.2*1000)/velocity[i]
                        y1 = ntrace
                    else:
                        raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
                    ax.plot([x0,x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5, label=str(velocity[i])+'m/s')
                    ax.plot([x0,-x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5)

                # set x axis
                ax.legend(loc='upper right', fontsize=7, shadow=False)
                ax.set_xlim(tt[0], tt[-1])
                ax.xaxis.set_ticks_position('bottom')
                ax.set_xlabel('Time(s)')

                # set y axis
                ax.invert_yaxis()
                yy_tick=np.linspace(0, ntrace, num=yticklabel_num)
                yy_label = np.around(np.linspace(dist_start, dist_end, num=yticklabel_num), decimals=2)
                ax.set_yticks(yy_tick)
                ax.set_yticklabels(yy_label)
                ax.set_ylabel('Distance(%s)' % (dist_unit))

                # set title
                if source_indx != None:
                    ax.set_title("StackData: source_indx=%d, win_num=%d, filter=[%.2f, %.2f] hz" % (source_indx, win_num, f1, f2))
                elif receiver_indx != None:
                    ax.set_title("StackData: receiver_indx=%d, win_num=%d, filter=[%.2f, %.2f] hz" % (source_indx, win_num, f1, f2))
                else:
                    raise ValueError("source_indx or receiver_indx must be specified")
                
                # save or show
                if save:
                    if save_path is None:
                        raise ValueError("save_path must be specified")
                    fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches='tight') 
                    plt.close(fig)
                else:
                    plt.show()
            else:
                raise ValueError("data is empty")
            
        else:
            raise ValueError("mode must be 'waveform' or 'mat'")



    def plot_moveout_all(self, 
                 corr_pair, 
                 pair_dist, 
                 dist_start=None, 
                 dist_end=None, 
                 amp_scale=1,
                 amp_normalize=True, 
                 win_num = 0, 
                 lag_start=None, 
                 lag_end=None, 
                 filter=False, 
                 f1=None, 
                 f2=None,
                 corners=4, 
                 zerophase=True,  
                 dist_interval=None, 
                 mode='waveform', 
                 cmap='seismic', 
                 linewidth=0.8,
                 yticklabel_num=10, 
                 figsize=(10, 6),
                 dist_unit="m", 
                 velocity=[], 
                 save=False, 
                 save_path=None, 
                 dpi=30):      
        
        self.maxlag = self.dt * (self.stack_data.shape[2]-1)/2
        
        # check corr_pair pair_dist
        if corr_pair.shape[0] != len(pair_dist):
            raise ValueError("the length of corr_pair is not consistent with pair_dist")
        
        # check demision of pair_dist
        if pair_dist.ndim == 2:
            pair_dist = pair_dist.reshape(-1)
        
        # check dist_unit
        if dist_unit not in ['m', 'km', 'degree']:
            raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
        
        # init dist_start
        if dist_start == None:
            dist_start = np.min(pair_dist)
        if dist_start < np.min(pair_dist):
            raise ValueError("dist_start must be >= (min(pair_dist)=%f)" % np.min(pair_dist))
        
        # init dist_end
        if dist_end == None:
            dist_end = np.max(pair_dist)
        if dist_end > np.max(pair_dist):
            raise ValueError("dist_end must be <= (max(pair_dist)=%f)" % np.max(pair_dist))

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

        # select dist_start dist_end
        id = np.where((pair_dist>=dist_start) & (pair_dist<=dist_end))[0]
        dist = pair_dist[id]

        # init dist_interval
        if dist_interval == None:
            dist_interval = np.min(dist)
        if dist_interval < np.diff(np.sort(dist)).min():
            raise ValueError("win_interval must be >= (np.diff(np.sort(dist)).min()=%f)" % np.diff(np.sort(dist)).min())
        if dist_interval > np.max(dist):
            raise ValueError("win_interval must be <= (np.max(dist)=%f)" % np.max(dist))
        
        # filter --> data
        if filter:
            data=np.empty((len(id), len(tt)))
            for i in range(0, len(id)):
                data[i,:] = bandpass(self.stack_data[id[i], win_num, tp_start:tp_end], f1, f2, int(1/self.dt), corners=corners, zerophase=zerophase)
        else:
            data = self.stack_data[id, win_num, tp_start:tp_end]
            
        # plot
        colors = list(mcolors.TABLEAU_COLORS.keys())
        if mode == 'waveform':
            if data.size != 0:
                # init fig
                fig, ax = plt.subplots(figsize=figsize)

                # normalize data
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
                
                # plot waveform
                for i in range(0, data_plot.shape[0]):
                    ax.plot(tt, amp_scale*data_plot[i]+dist[i], 'k', linewidth=linewidth)

                # plot velocity
                for i in range(0, len(velocity)):
                    x0=0; y0=0
                    if dist_unit == 'm':
                        x1 = (dist_end-dist_start)/velocity[i]
                        y1 = dist_end
                    elif dist_unit == 'km':
                        x1 = (dist_end-dist_start)/1000/velocity[i]
                        y1 = dist_end
                    elif dist_unit == 'degree':
                        x1 = (dist_end-dist_start)/(111.2*1000)/velocity[i]
                        y1 = dist_end
                    else:
                        raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
                    ax.plot([x0,x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5, label=str(velocity[i])+'m/s')
                    ax.plot([x0,-x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5)

                # set axis
                ax.set_xlim(tt[0], tt[-1])
                ax.set_xlabel('Time(s)')
                ax.set_ylabel('Distance(%s)' % (dist_unit))
                ax.legend(loc='upper right', fontsize=8, shadow=False)
                ax.set_title("StackData: win_num=%d, filter=[%.2f, %.2f] hz" % (win_num, f1, f2))
            
                # save or show
                if save:
                    if save_path is None:
                        raise ValueError("save_path must be specified")
                    fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches = 'tight') 
                    plt.close(fig)
                else:
                    plt.show()

            else:
                raise ValueError("data is empty")

        elif mode=="mat":
            if data.size != 0:
                # init fig
                fig, ax = plt.subplots(figsize=figsize)

                # average waveforms
                ntrace = int((dist_end-dist_start)/dist_interval)
                ndata  = np.zeros(shape=(ntrace,len(tt)))
                for i in range(0,ntrace):
                    tindx = np.where(((dist-dist_start)>=i*dist_interval) & ((dist-dist_start)<(i+1)*dist_interval))[0]
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
                
                # plot
                cax=ax.matshow(ndata_plot, extent=[lag_start, lag_end, ntrace, 0], aspect='auto', cmap=cmap) # cmap=cmaps.MPL_jet,
                for i in range(0, len(velocity)):
                    x0=0
                    y0=0
                    if dist_unit == 'm':
                        x1 = (dist_end-dist_start)/velocity[i]
                        y1 = ntrace
                    elif dist_unit == 'km':
                        x1 = (dist_end-dist_start)/1000/velocity[i]
                        y1 = ntrace
                    elif dist_unit == 'degree':
                        x1 = (dist_end-dist_start)/(111.2*1000)/velocity[i]
                        y1 = ntrace
                    else:
                        raise ValueError("dist_unit must be 'm', 'km', or 'degree'")
                    ax.plot([x0,x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5, label=str(velocity[i])+'m/s')
                    ax.plot([x0,-x1], [y0,y1], color=colors[i], linestyle='--', linewidth=1.5)

                # set x axis
                ax.set_title("StackData: win_num=%d, filter=[%.2f, %.2f] hz" % (win_num, f1, f2))
                ax.legend(loc='upper right', fontsize=7, shadow=False)
                ax.set_xlim(tt[0], tt[-1])
                ax.xaxis.set_ticks_position('bottom')
                ax.set_xlabel('Time(s)')

                # set y axis
                ax.invert_yaxis()
                yy_tick=np.linspace(0, ntrace, num=yticklabel_num)
                yy_label = np.around(np.linspace(dist_start, dist_end, num=yticklabel_num), decimals=2)
                ax.set_yticks(yy_tick)
                ax.set_yticklabels(yy_label)
                ax.set_ylabel('Distance(%s)' % (dist_unit))

                # save or show
                if save:
                    if save_path is None:
                        raise ValueError("save_path must be specified")
                    fig.savefig(os.path.join(save_path), dpi=dpi, format='pdf', bbox_inches='tight') 
                    plt.close(fig)
                else:
                    plt.show()
            else:
                raise ValueError("data is empty")
            
        else:
            raise ValueError("mode must be 'waveform' or 'mat'")


