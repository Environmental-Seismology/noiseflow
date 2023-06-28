/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef FILTER_HPP
#define FILTER_HPP

#include "utils.hpp"

namespace SIGNAL {

using namespace kfr;

// bandpass_filter
template <class T1>
std::unique_ptr<biquad_filter<T1>> bandpass_filter(double freqmin, double freqmax, double df, int corners) {
    double fe = 0.5 * df;
    double low = freqmin / fe;
    double high = freqmax / fe;

    std::unique_ptr<biquad_filter<T1>> filter;
    std::vector<biquad_params<T1>> bqs = to_sos(iir_bandpass(butterworth<T1>(corners), low, high));
    filter.reset(new biquad_filter<T1>(bqs));

    return filter;
}


template <class T1, class T2>
void bandpass(T2 &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads) {
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    std::vector<std::size_t> shape = {std::size_t(samples_num)};

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start bandpass with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) { 
    std::cout << "Start bandpass with " << 1 << " threads without openmp..." << std::endl;
  }
#endif
    for(int i=0; i<channels_num; i++){
    // flag
#ifdef _OPENMP
    if (flag && omp_get_thread_num()==0 && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (flag && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif
        std::unique_ptr<biquad_filter<T1>> filter = bandpass_filter<T1>(freqmin, freqmax, df, corners);
        auto input = make_univector(data.data()+i*samples_num, samples_num);
        univector<T1> output(samples_num);
        if (zerophase) {
            univector<T1> output_first(samples_num);
            filter->apply(output_first, input);
            std::reverse(output_first.begin(), output_first.end());
            filter->apply(output, make_univector(output_first));
            std::reverse(output.begin(), output.end());
        }
        else {
            filter->apply(output, input);
        }
        xt::view(data, i, xt::all()) = xt::adapt(output.data(), samples_num, xt::no_ownership(), shape);
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End bandpass with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}



// lowpass_filter
template <class T1>
std::unique_ptr<biquad_filter<T1>> lowpass_filter(double freqlow, double df, int corners) {
    std::unique_ptr<biquad_filter<T1>> filter;
    std::vector<biquad_params<T1>> bqs = to_sos(iir_lowpass(butterworth<T1>(corners), freqlow, df));
    filter.reset(new biquad_filter<T1>(bqs));

    return filter;
}


template <class T1, class T2>
void lowpass(T2 &data, double freqlow, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads) {
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    std::vector<std::size_t> shape = {std::size_t(samples_num)};

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start lowpass with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start lowpass with " << 1 << " threads without openmp..." << std::endl;
  }
#endif
    for(int i=0; i<channels_num; i++){
    // flag
#ifdef _OPENMP
    if (flag && omp_get_thread_num()==0 && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (flag && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif
        std::unique_ptr<biquad_filter<T1>> filter = lowpass_filter<T1>(freqlow, df, corners);
        auto input = make_univector(data.data()+i*samples_num, samples_num);
        univector<T1> output(samples_num);
        if (zerophase) {
            univector<T1> output_first(samples_num);
            filter->apply(output_first, input);
            std::reverse(output_first.begin(), output_first.end());
            filter->apply(output, make_univector(output_first));
            std::reverse(output.begin(), output.end());
        }
        else {
            filter->apply(output, input);
        }
        xt::view(data, i, xt::all()) = xt::adapt(output.data(), samples_num, xt::no_ownership(), shape);
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End lowpass with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}



// highpass_filter
template <class T1>
std::unique_ptr<biquad_filter<T1>> highpass_filter(double freqhigh, double df, int corners) {
    std::unique_ptr<biquad_filter<T1>> filter;
    std::vector<biquad_params<T1>> bqs = to_sos(iir_highpass(butterworth<T1>(corners), freqhigh, df));
    filter.reset(new biquad_filter<T1>(bqs));

    return filter;
}


template <class T1, class T2>
void highpass(T2 &data, double freqhigh, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads) {
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    std::vector<std::size_t> shape = {std::size_t(samples_num)};

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start highpass with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start highpass with " << 1 << " threads without openmp..." << std::endl;
  }
#endif
    for(int i=0; i<channels_num; i++){
    // flag
#ifdef _OPENMP
    if (flag && omp_get_thread_num()==0 && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (flag && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif
        std::unique_ptr<biquad_filter<T1>> filter = highpass_filter<T1>(freqhigh, df, corners);
        auto input = make_univector(data.data()+i*samples_num, samples_num);
        univector<T1> output(samples_num);
        if (zerophase) {
            univector<T1> output_first(samples_num);
            filter->apply(output_first, input);
            std::reverse(output_first.begin(), output_first.end());
            filter->apply(output, make_univector(output_first));
            std::reverse(output.begin(), output.end());
        }
        else {
            filter->apply(output, input);
        }
        xt::view(data, i, xt::all()) = xt::adapt(output.data(), samples_num, xt::no_ownership(), shape);
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End highpass with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}




// bandstop_filter
template <class T1>
std::unique_ptr<biquad_filter<T1>> bandstop_filter(double freqmin, double freqmax, double df, int corners) {
    double fe = 0.5 * df;
    double low = freqmin / fe;
    double high = freqmax / fe;

    std::unique_ptr<biquad_filter<T1>> filter;
    std::vector<biquad_params<T1>> bqs = to_sos(iir_bandstop(butterworth<T1>(corners), low, high));
    filter.reset(new biquad_filter<T1>(bqs));

    return filter;
}


template <class T1, class T2>
void bandstop(T2 &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads) {
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    std::vector<std::size_t> shape = {std::size_t(samples_num)};

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start bandstop with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start bandstop with " << 1 << " threads without openmp..." << std::endl;
  }
#endif
    for(int i=0; i<channels_num; i++){
    // flag
#ifdef _OPENMP
    if (flag && omp_get_thread_num()==0 && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (flag && i%flag_gap==0 && flag_gap!=channels_num) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif
        std::unique_ptr<biquad_filter<T1>> filter = bandstop_filter<T1>(freqmin, freqmax, df, corners);
        auto input = make_univector(data.data()+i*samples_num, samples_num);
        univector<T1> output(samples_num);
        if (zerophase) {
            univector<T1> output_first(samples_num);
            filter->apply(output_first, input);
            std::reverse(output_first.begin(), output_first.end());
            filter->apply(output, make_univector(output_first));
            std::reverse(output.begin(), output.end());
        }
        else {
            filter->apply(output, input);
        }
        xt::view(data, i, xt::all()) = xt::adapt(output.data(), samples_num, xt::no_ownership(), shape);
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End bandstop with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}

// **************************************************************************
// *                        end
// **************************************************************************

}

#endif 