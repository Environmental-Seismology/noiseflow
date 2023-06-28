import numpy as np

from obspy.core import UTCDateTime


# get rfft/corr timestamp vector
def get_timestamp(win_num,             
                cc_len, 
                cc_step,
                t_min=UTCDateTime("1970-01-01T00:00:00.0")):

    t_max = t_min + (win_num-1)*(cc_len-cc_step) + cc_len
    win_interval = cc_len-cc_step
    win_start_vector = np.array([t_min + i*win_interval for i in range(0, win_num)])

    return t_max, win_interval, win_start_vector


# get stack timestamp vector
def get_stack_timestamp(win_num,
                        stack_len,
                        stack_step,
                        cc_len, 
                        cc_step,
                        t_min=UTCDateTime("1970-01-01T00:00:00.0")):
    
    corr_win_num = (win_num-1)*(stack_len-stack_step) + stack_len
    t_max = t_min + (corr_win_num-1)*(cc_len-cc_step) + cc_len
    win_len = stack_len*(cc_len-cc_step)
    win_interval = (stack_len-stack_step)*(cc_len-cc_step)
    win_start_vector = np.array([t_min + i*win_interval for i in range(0, win_num)])

    return t_max, win_len, win_interval, win_start_vector
    

def time_linspace(start_time, end_time, num):
    t_min = 0
    t_max = end_time - start_time
    t_vector = np.linspace(t_min, t_max, num=num)
    date_vector = np.array([start_time + i for i in t_vector])
              
    return date_vector
