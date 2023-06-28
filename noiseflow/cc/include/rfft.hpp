/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef RFFT_HPP
#define RFFT_HPP

#include "utils.hpp"

namespace CC {

// Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
template <typename T1, typename T2, typename T3> 
class RFFTClass{
  private:
    T2 _raw_data; // input data
    float _dt;
    float _cc_len;
    float _cc_step;
    std::string _time_norm;
    float _clip_std;
    int _smooth_N;
    std::string _freq_norm; 
    float _freqmin;
    float _freqmax;
    int _whiten_npad;
    int _smoothspect_N;
    bool _flag;
    int _flag_gap;
    int _threads;
    
  public:
    int _cc_len_points;
    int _cc_step_points;
    int _npts; // number of points for each channel in the rawdata 
    int _rfft_npts; // number of points for each channel in the rfftdata
    int _channel_num; // number of channels in the rawdata 
    int _win_num; // number of windows for each channel in the rawdata
    xt::xtensor<int, 2> _win_info; // slice window
    T3 _output_data; // output data

    // constructor
    RFFTClass(T2& raw_data, float dt, float cc_len, float cc_step, std::string time_norm, float clip_std, int smooth_N, std::string freq_norm, float freqmin, float freqmax, int whiten_npad, int smoothspect_N, bool flag, int flag_gap, int threads);

    // copy constructor
    // RFFTClass(RFFTClass & A);

    // destructor
    ~RFFTClass();

    // time_norm --> time_norm_data 
    xt::xarray<T1>  time_norm_no(int channel);
    xt::xarray<T1>  time_norm_onebit(int channel);
    xt::xarray<T1>  time_norm_clip(int channel);
    xt::xarray<T1>  time_norm_smooth(int channel);

    // rfft --> rfft_data
    xt::xarray<std::complex<T1>> rfft(xt::xarray<T1>& time_norm_data);

    // freq_norm --> rfft_norm_data
    xt::xarray<std::complex<T1>> freq_norm_no(xt::xarray<std::complex<T1>>& rfft_data);
    xt::xarray<std::complex<T1>> freq_norm_whiten(xt::xarray<std::complex<T1>>& rfft_data);

    // run workflow
    void run();

    // return
    T3& get_rfftdata() {return _output_data;};
};




// **************************************************************************
// *                        constructor
// **************************************************************************
// parameter constructor                                                                                                                                                      
template <typename T1, typename T2, typename T3> 
RFFTClass<T1, T2, T3>::RFFTClass(T2& raw_data, float dt, float cc_len, float cc_step, std::string time_norm, float clip_std, int smooth_N, std::string freq_norm, float freqmin, float freqmax, int whiten_npad, int smoothspect_N, bool flag, int flag_gap, int threads) 
{
  _raw_data = raw_data;
  _dt = dt;
  _cc_len = cc_len;
  _cc_step = cc_step;
  _time_norm = time_norm;
  _clip_std = clip_std;
  _smooth_N = smooth_N;
  _freq_norm = freq_norm;
  _freqmin = freqmin;
  _freqmax = freqmax;
  _whiten_npad = whiten_npad;
  _smoothspect_N = smoothspect_N;
  _flag = flag;
  _flag_gap = flag_gap;
  _threads = threads;

  // convert time to points
  _cc_len_points = int(_cc_len/_dt);
  _cc_step_points = int(_cc_step/_dt);
  _npts = raw_data.shape(1);
  _rfft_npts = _cc_len_points/2+1;
  _channel_num = raw_data.shape(0);
  
  // slice_window --> _win_num, _win_info
  _win_info = CC::slice_window(_npts, _cc_len_points, _cc_step_points);
  _win_num = _win_info.shape(0);
}

                                                                                                                                                       
// virtual destructor                                                                                                                                                       
template <typename T1, typename T2, typename T3> 
RFFTClass<T1, T2, T3>::~RFFTClass() {}




// **************************************************************************
// *                        time_norm
// **************************************************************************
// no
template <typename T1, typename T2, typename T3> 
inline xt::xarray<T1> RFFTClass<T1, T2, T3>::time_norm_no(int channel)
{
  xt::xarray<T1> time_norm_data = xt::view(_raw_data, channel, xt::all());

  return time_norm_data;
}


// onebit
template <typename T1, typename T2, typename T3> 
inline xt::xarray<T1> RFFTClass<T1, T2, T3>::time_norm_onebit(int channel)
{
  xt::xarray<T1> time_norm_data = xt::sign(xt::view(_raw_data, channel, xt::all()));
  
  return time_norm_data;
}


// clip
template <typename T1, typename T2, typename T3> 
inline xt::xarray<T1> RFFTClass<T1, T2, T3>::time_norm_clip(int channel)
{
  xt::xarray<T1> time_norm_data = xt::view(_raw_data, channel, xt::all());
  T1 lim = _clip_std * xt::stddev(time_norm_data)(0);
  filtration(time_norm_data, time_norm_data>lim) = lim;
  filtration(time_norm_data, time_norm_data<(-lim)) = -lim;

  return time_norm_data;
}


// smooth
template <typename T1, typename T2, typename T3> 
inline xt::xarray<T1> RFFTClass<T1, T2, T3>::time_norm_smooth(int channel)
{
  xt::xarray<T1> data = xt::view(_raw_data, channel, xt::all());
  xt::xarray<T1> time_norm_data = data / CC::moving_ave<T1>(xt::abs(data), _smooth_N);

  return time_norm_data;
}





// **************************************************************************
// *                        rfft
// **************************************************************************
template <typename T1, typename T2, typename T3> 
inline xt::xarray<std::complex<T1>> RFFTClass<T1, T2, T3>::rfft(xt::xarray<T1>& time_norm_data)
{
  xt::xarray<std::complex<T1>> rfft_data= xt::empty<T1>({_win_num, _rfft_npts});

  for(int i=0; i<_win_num; i++){
    auto segment = xt::view(time_norm_data, xt::range(_win_info(i,0), _win_info(i,1)), xt::all());
    auto segment_rfft = xt::fftw::rfft(xt::eval(segment));
    xt::view(rfft_data, i, xt::all()) = segment_rfft; 
  }

  return rfft_data;
}





// **************************************************************************
// *                        freq_norm
// **************************************************************************
template <typename T1, typename T2, typename T3> 
inline xt::xarray<std::complex<T1>> RFFTClass<T1, T2, T3>::freq_norm_no(xt::xarray<std::complex<T1>>& rfft_data)
{
  auto freq_norm_data = rfft_data;

  return freq_norm_data;
}


template <typename T1, typename T2, typename T3> 
inline xt::xarray<std::complex<T1>> RFFTClass<T1, T2, T3>::freq_norm_whiten(xt::xarray<std::complex<T1>>& rfft_data)
{
  xt::xarray<std::complex<T1>> rfft_norm_data = CC::whiten<T1>(rfft_data, _dt, _freq_norm, _freqmin, _freqmax, _smoothspect_N, _whiten_npad);
  
  return rfft_norm_data;
}






// **************************************************************************
// *                        run
// **************************************************************************
template <typename T1, typename T2, typename T3> 
void RFFTClass<T1, T2, T3>::run()
{
  auto start_time = std::chrono::high_resolution_clock::now();
  auto end_time = start_time;
  _output_data = xt::empty<std::complex<T1>>({_channel_num, _win_num, _rfft_npts});
  
#ifdef _OPENMP
  omp_set_num_threads(_threads);
  std::cout << "Start rfft with " << _threads  << " threads using openmp..." << std::endl;
  #pragma omp parallel for
#else
  std::cout << "Start rfft with " << 1 << " threads without openmp..." << std::endl;
#endif

  for(int i = 0; i < _channel_num; i++){
// flag
#ifdef _OPENMP
    if (_flag && omp_get_thread_num()==0 && i%_flag_gap==0) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (_flag && i%_flag_gap==0) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process channel " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif

    // time_norm
    xt::xarray<T1> time_norm_data;
    if (_time_norm == "no") {
      time_norm_data = time_norm_no(i);
    } 
    else if (_time_norm == "clip"){
      time_norm_data = time_norm_clip(i);
    } 
    else if (_time_norm == "onebit") {
      time_norm_data = time_norm_onebit(i);
    } 
    else if (_time_norm == "smooth") {
      time_norm_data = time_norm_smooth(i);
    } 
    else {
      std::cout << "error: please input the correct time_norm method." << std::endl;
      exit(1);
    }

    // rfft
    xt::xarray<std::complex<T1>> rfft_data = rfft(time_norm_data);

    // freq_norm
    xt::xarray<std::complex<T1>> rfft_norm_data;
    if (_freq_norm == "no") {
      rfft_norm_data = freq_norm_no(rfft_data);
    } else {
      rfft_norm_data = freq_norm_whiten(rfft_data);
    }

    // fill _output_data
    xt::view(_output_data, i, xt::all()) = rfft_norm_data; 
  }

  end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
  std::cout << "End rfft with total time " << elapsed_time.count() << "s" <<std::endl;

}



// **************************************************************************
// *                        end
// **************************************************************************

}

#endif 