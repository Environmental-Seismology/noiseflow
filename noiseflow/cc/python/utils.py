import numpy as np


def split(num, n_jobs):
    chunk_size = num // n_jobs
    results = []
    for i in range(0, n_jobs):
        chunk_start = i*chunk_size
        if i == n_jobs - 1:
            chunk_end = num
        else:
            chunk_end = (i+1)*chunk_size
        results.append([chunk_start, chunk_end])

    return results


def slice_window(npts, segment_points, step_points):
    if segment_points < npts:
        slide_points = segment_points - step_points
        win_num = 0
        for i in range(0, int(npts/slide_points)):
            if (i * slide_points + segment_points) <= npts:
                win_num += 1
            else:
                break
        win_info = np.empty((win_num, 2), dtype=int)
        for i in range(win_num):
            win_info[i, 0] = i * slide_points
            win_info[i, 1] = i * slide_points + segment_points
    elif segment_points == npts:
        win_num = 1
        win_info = np.array([[0, npts]], dtype=int)
    else:
        raise ValueError("error: segment-points length is larger than npts when slicing windows!")
    
    return win_info


def moving_ave(A, N):
    temp = np.zeros(A.shape[0] + 2*N)

    temp_len = temp.shape[0]
    temp[N:temp_len-N] = A
    temp[0:N] = temp[N]
    temp[temp_len-N:temp_len] = temp[temp_len-N-1]

    nn = np.ones(N) / N
    b1 = np.convolve(temp, nn, mode='full')

    n1 = N + (N-1) // 2
    n2 = N + (N-1 - (N-1) // 2)
    B1 = b1[n1 : b1.shape[0] - n2]

    return B1


def whiten(rfft_data, dt, freq_norm_method, freqmin, freqmax, smoothspect_N, whiten_npad):
    _i = 0.0 + 1.0j
    _i0 = 0.0 + 0.0j
    _pi = np.pi

    win_num = rfft_data.shape[0]
    fft_npts = rfft_data.shape[1]

    freq_array = np.arange(fft_npts) / dt / 2.0 / (fft_npts-1) # note: must divide by 2.0, because rfft_npts is not npts.
    J = np.where((freq_array>=freqmin) & (freq_array<=freqmax))[0]

    low = J[0] - whiten_npad
    if (low <= 0):
        low = 0
    left = J[0]
    right = J[-1]
    high = J[-1] + whiten_npad
    if (high > fft_npts):
        high = fft_npts

    rfft_norm_data = _i0*np.ones((win_num, fft_npts))

    for i in range(0, win_num):
        #  left zero cut-off
        # // xt::view(rfft_norm_data, i, xt::range(0, low)) = _i0*xt::view(rfft_data, i, xt::range(0, low));

        # // left tapering
        smo1 = np.power(np.cos( np.linspace(_pi/2, _pi, left-low)), 2)
        exp1 = np.exp(_i*np.angle(rfft_data[i, low:left]))
        rfft_norm_data[i, low:left] = smo1*exp1
 
        # // pass band
        if freq_norm_method == "whiten":
            smo2 = np.ones((right-left))
            exp2 = np.exp(_i*np.angle(rfft_data[i, left: right]))
            rfft_norm_data[ i, left:right] = smo2*exp2
        elif freq_norm_method == "smooth_whiten":
            data = rfft_data[i, left:right]
            rfft_norm_data[ i, left:right] = data / moving_ave(np.abs(data), smoothspect_N)
    
        # // right tapering
        smo3 = np.power(np.cos( np.linspace(0, _pi/2, high-right)), 2)
        exp3 = np.exp(_i*np.angle(rfft_data[ i, right: high]))
        rfft_norm_data[i, right: high] = smo3*exp3

        # // right zero cut-off
        # // xt::view(rfft_norm_data, i, xt::range(high, -1)) = _i0*xt::view(rfft_data, i, xt::range(high, -1));
    
    return rfft_norm_data

