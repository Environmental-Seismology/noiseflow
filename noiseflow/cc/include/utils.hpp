/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef UTILS_HPP
#define UTILS_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <vector>
#include <string>
// #include <ctime>
#include <chrono>
#include <complex>
#include <exception>

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xmanipulation.hpp> // for xt::roll
#include <xtensor/xadapt.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xindex_view.hpp>
#include <xsimd/xsimd.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor-fftw/basic.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>


namespace CC {

inline xt::xtensor<int, 2> slice_window(int npts, int segment_points, int step_points)
{
    int win_num;
    xt::xtensor<int, 2> win_info;

    if (segment_points < npts){
        int slide_points = segment_points - step_points;
        win_num = 0;
        for (int i=0; i<(npts/slide_points); i++){
            if ((i * slide_points + segment_points) <= npts){
                win_num++;
            } 
            else {
                break;
            }
        }
        win_info = xt::empty<int>({win_num, 2});
        for (int i=0; i<win_num; i++){
            win_info(i,0) = i * slide_points;
            win_info(i,1) = i * slide_points + segment_points;
        }
    }
    else if (segment_points == npts){
        win_num = 1;
        win_info = xt::xtensor<int, 2> {{0, npts}};
    } 
    else if (segment_points > npts){
        std::cout << "error: segment-points length is larger than npts when slicing windows!" << std::endl;
        exit(1);
    } 
    else {
        std::cout << "error: segment-points length is small than 0!" << std::endl; 
        exit(1);
    }

    return win_info;
}

  

// A: 1D xarray
template <class T>
inline xt::xarray<T> moving_ave(const xt::xarray<T>& A, int N)
{
    xt::xarray<T> temp = xt::zeros<T>({A.shape(0)+2*N});
    
    int temp_len = temp.shape(0);
    xt::view(temp, xt::range(N, temp_len-N)) = A;
    xt::view(temp, xt::range(0, N)) = temp(N);
    xt::view(temp, xt::range(temp_len-N, temp_len)) = temp(temp_len-N-1);

    xt::xarray<T> nn = xt::ones<T>({N}) / T(N);
    auto b1 = xt::convolve(temp, nn, xt::convolve_mode::full());

    int n1 = N+(N-1)/2;
    int n2 = N+(N-1-(N-1)/2);
    auto B1 = xt::view(b1, xt::range(n1, b1.shape(0)-n2));

    return B1;
}


// rfft_data: 2D xarray, shape: (win_num, rfft_npts)
template <class T>
inline xt::xarray<std::complex<T>> whiten(const xt::xarray<std::complex<T>>& rfft_data, float dt, std::string freq_norm, float freqmin, float freqmax, int smoothspect_N, int whiten_npad)
{
    std::complex<T> _i {0, 1};
    std::complex<T> _i0 {0, 0};
    T _pi = xt::numeric_constants<T>::PI;

    int win_num = rfft_data.shape(0);
    int rfft_npts = rfft_data.shape(1);

    xt::xarray<T> freq_array = xt::arange(rfft_npts) / dt / 2.0 / T(rfft_npts-1); // note: must divide by 2.0, because rfft_npts is not _npts.
    auto J = xt::where((freq_array>=freqmin) & (freq_array<=freqmax))[0];

    int low = J[0] - whiten_npad;
    if (low <= 0)
        low = 0;
    int left = J[0];
    int right = J[J.size()-1];
    int high = J[J.size()-1] + whiten_npad;
    if (high > rfft_npts)
        high = rfft_npts;

    xt::xarray<std::complex<T>> rfft_norm_data = _i0*xt::ones<std::complex<T>>({win_num, rfft_npts});

    for(int i=0; i<win_num; i++){
        // left zero cut-off
        // xt::view(rfft_norm_data, i, xt::range(0, low)) = _i0*xt::view(rfft_data, i, xt::range(0, low));

        // left tapering
        auto smo1 = xt::pow(xt::cos( xt::linspace<T>(_pi/2, _pi, left-low)), 2);
        auto exp1 = xt::exp(_i*xt::arg(xt::view(rfft_data, i, xt::range(low, left))));
        xt::view(rfft_norm_data, i, xt::range(low, left)) = smo1*exp1;
 
        // pass band
        if (freq_norm == "whiten"){
            auto smo2 = xt::ones<T>({right-left});
            auto exp2 = xt::exp(_i*xt::arg(xt::view(rfft_data, i, xt::range(left, right))));
            xt::view(rfft_norm_data, i, xt::range(left, right)) = smo2*exp2;
        } 
        else if(freq_norm == "smooth_whiten"){
            auto data = xt::view(rfft_data, i, xt::range(left, right));
            xt::view(rfft_norm_data, i, xt::range(left, right)) = data / moving_ave<std::complex<T>>(xt::abs(data), smoothspect_N);
        } 
        else {
            std::cout << "freq_norm is not defined" << std::endl;
            exit(1);
        }
    
        // right tapering
        auto smo3 = xt::pow(xt::cos( xt::linspace<T>(0, _pi/2, high-right)), 2);
        auto exp3 = xt::exp(_i*xt::arg(xt::view(rfft_data, i, xt::range(right, high))));
        xt::view(rfft_norm_data, i, xt::range(right, high)) = smo3*exp3;

        // right zero cut-off
        // xt::view(rfft_norm_data, i, xt::range(high, -1)) = _i0*xt::view(rfft_data, i, xt::range(high, -1));
    }

    return rfft_norm_data;
}



// data: (nchannel, npts) 2D array
template <class T>
xt::xarray<T> robust(const xt::xarray<T>& data)
{
    int nstep = 0;
    int maxstep = 10;
    int channel_num = data.shape(0);

    T res = 9e9;
    T epsilon = 1e-5;
    
    xt::xarray<T> w = xt::ones<T>({channel_num});
    xt::xarray<T> newstack = xt::median(data,0);
    xt::xarray<T> crap;
    xt::xarray<T> crap_dot;
    xt::xarray<T> di_norm;
    xt::xarray<T> ri;
    xt::xarray<T> ri_norm;
    xt::xarray<T> stack;

    while (res>epsilon && nstep <=maxstep){
        stack = newstack;

        for(int i=0; i<channel_num;i++){
            crap = stack * xt::transpose(xt::view(data, i, xt::all()));
            crap_dot = xt::sum(crap);
            di_norm = xt::linalg::norm(xt::view(data, i, xt::all()));
            ri = xt::view(data, i, xt::all()) - crap_dot*stack;
            ri_norm = xt::linalg::norm(ri);
            xt::view(w, i) = xt::abs(crap_dot)/di_norm/ri_norm;
        }

        w = w/xt::sum(w);
        newstack = xt::sum(xt::transpose(w*xt::transpose(data)), 0);
        res = xt::linalg::norm(newstack-stack, 1)/xt::linalg::norm(newstack)/T(channel_num);
        nstep += 1;
    }

    return newstack;
}



// x: (npts) 1D complex-array
template <class T>
inline xt::xarray<std::complex<T>> hilbert_transform(const xt::xarray<std::complex<T>>& x)
{
    int N = x.shape(0);
    xt::xarray<std::complex<T>> Xf = xt::fftw::fft(x);
    xt::xarray<T> h = xt::zeros<T>({N});

    if ((N%2) == 0){
        h(0) = 1.0;
        h(N/2) = 1.0;
        xt::view(h, xt::range(1,N/2)) = 2.0;
    } 
    else{
        h(0) = 1.0;
        xt::view(h, xt::range(1,(N+1)/2)) = 2.0;
    }
        
    return xt::fftw::ifft(xt::eval(Xf * h));
}



// x: (nchannel, npts) 2D array
template <class T>
inline xt::xarray<T> pws(const xt::xarray<T>& x)
{   
    int channel_num = x.shape(0);

    if (channel_num == 0){
        xt::xarray<T> newstack;
        return newstack;
    } 
    else{
        int row_num = x.shape(0);
        int colume_num = x.shape(1);
        std::complex<float> _i {0, 1};
        xt::xarray<std::complex<T>> analytic = xt::empty<std::complex<T>>({row_num, colume_num});

        for (int i=0; i<row_num; i++){
            xt::view(analytic, i, xt::all()) = hilbert_transform<T>(xt::view(x, i, xt::all())); 
        }

        auto phase = xt::angle(analytic);
        auto phase_stack = xt::mean(xt::exp(_i*phase), 0); 
        // auto phase_stack = xt::sum(xt::exp(_i*phase),0) / static_cast<T>(row_num);

        auto phase_stack_new = xt::pow(xt::abs(phase_stack), 2);
        xt::xarray<T> weighted = xt::mean(x*phase_stack_new, 0);
        // xt::xarray<T> weighted = xt::sum(x*phase_stack_new, 0) / static_cast<T>(row_num);

        return weighted;
    }
}


} 
#endif