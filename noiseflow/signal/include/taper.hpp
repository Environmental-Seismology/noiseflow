/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef TAPER_HPP
#define TAPER_HPP

#include "utils.hpp"

namespace SIGNAL {

using namespace kfr;



// taper
template <class T1, class T2>
void taper(T2 &data, double max_percentage, std::string type, std::string side, bool flag, int flag_gap, int threads) 
{
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    int wlen = std::floor(samples_num * max_percentage);

    univector<T1> windows;
    if (type == "hann"){
        windows = window_hann(wlen*2+1);
    }
    else if (type == "hamming"){
        windows = window_hamming(wlen*2+1);
    }
    else if (type == "blackman"){
        windows = window_blackman(wlen*2+1);
    }
    else if (type == "blackman_harris"){
        windows = window_blackman_harris(wlen*2+1);
    }
    else if (type == "gaussian"){
        windows = window_gaussian(wlen*2+1);
    }
    else if (type == "triangular"){
        windows = window_triangular(wlen*2+1);
    }
    else if (type == "bartlett"){
        windows = window_bartlett(wlen*2+1);
    }
    else if (type == "cosine"){
        windows = window_cosine(wlen*2+1);
    }
    else if (type == "cosine_np"){
        windows = window_cosine_np(wlen*2+1);
    }
    else if (type == "bartlett_hann"){
        windows = window_bartlett_hann(wlen*2+1);
    }
    else if (type == "bohman"){
        windows = window_bohman(wlen*2+1);
    }
    else if (type == "lanczos"){
        windows = window_lanczos(wlen*2+1);
    }
    else if (type == "flattop"){
        windows = window_flattop(wlen*2+1);
    }
    else if (type == "window_kaiser"){
        windows = window_kaiser(wlen*2+1, 2.5);
    }
    else {
        std::cout << "error: please input the correct windows method." << std::endl;
    }

    std::vector<std::size_t> shape = {std::size_t(wlen*2+1)};
    auto taper_sides = xt::adapt(windows.data(), wlen*2+1, xt::no_ownership(), shape);

    xt::xarray<T1> taper;
    if (side=="left"){
        taper = xt::hstack(xtuple(xt::view(taper_sides, xt::range(0, wlen)), xt::ones<T1>({samples_num-wlen})));
    }
    else if (side=="right"){
        taper = xt::hstack(xtuple(xt::ones<T1>({samples_num-wlen}), xt::view(taper_sides, xt::range(wlen+1, wlen*2+1))));
    }
    else if (side=="both"){
        taper = xt::hstack(xtuple(xt::view(taper_sides, xt::range(0, wlen)), xt::ones<T1>({samples_num-wlen*2}), xt::view(taper_sides, xt::range(wlen+1, wlen*2+1))));
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start taper with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start taper with " << 1 << " threads without openmp..." << std::endl;
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
        xt::view(data, i, xt::all()) *= taper;
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End taper with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}


// **************************************************************************
// *                        end
// **************************************************************************

}

#endif 