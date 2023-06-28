/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef DETREND_HPP
#define DETREND_HPP

#include "utils.hpp"

namespace SIGNAL {

using namespace kfr;



// detrend
template <class T1, class T2>
void detrend(T2 &data, std::string type, bool flag, int flag_gap, int threads) 
{
    int samples_num = data.shape(1);
    int channels_num = data.shape(0);
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;

#ifdef _OPENMP
  omp_set_num_threads(threads);
  if (flag) {
    std::cout << "Start detrend with " << threads  << " threads using openmp..." << std::endl;
  }
  #pragma omp parallel for
#else
  if (flag) {
    std::cout << "Start detrend with " << 1 << " threads without openmp..." << std::endl;
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
        if (type == "linear"){
            xt::xtensor<T1, 2> A = xt::ones<T1>({samples_num, 2});
            xt::view(A, xt::all(), 0) = xt::arange(samples_num);
            auto y = xt::view(data, i, xt::all());
            auto result = xt::linalg::lstsq(A, y);
            auto trend = xt::linalg::dot(A, std::get<0>(result));
            xt::view(data, i, xt::all()) -= trend;
        }
        else if (type == "demean"){
            auto dd = xt::view(data, i, xt::all());
            auto trend = xt::mean(dd);
            xt::view(data, i, xt::all()) -= trend;
        }
        else {
            std::cout << "error: please input the correct detrend method." << std::endl;
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if (flag) {
        std::cout << "End detrend with total time " << elapsed_time.count() << "s" <<std::endl;
    }
}



// **************************************************************************
// *                        end
// **************************************************************************

}

#endif 