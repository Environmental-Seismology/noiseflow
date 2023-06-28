/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef DECIMATE_HPP
#define DECIMATE_HPP

#include "utils.hpp"

namespace SIGNAL {

using namespace kfr;



// decimate
template <class T1, class T2>
T2 decimate(T2 &data, double df, double df_new, bool flag, int flag_gap, int threads) 
{
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    int decimates_num = std::floor(samples_num * df_new/df);
    int interval = std::round(df/df_new);
    T2 re_data= xt::empty<T1>({channels_num, decimates_num});
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start decimate with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start decimate with " << 1 << " threads without openmp..." << std::endl;
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

        for (int j=0; j<decimates_num; j++){
            re_data(i, j) = data(i, j*interval);
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End decimate with total time " << elapsed_time.count() << "s" <<std::endl;
    }
    return re_data;
}

// **************************************************************************
// *                        end
// **************************************************************************

}

#endif 