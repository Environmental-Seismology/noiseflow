/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef CORR_HPP
#define CORR_HPP

#include "utils.hpp"

namespace CC {

// Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
template <typename T1, typename T2, typename T3> 
class CorrClass{
  private:
    T2 _fft_data; // input data
    float _dt;
    std::string _corr_method;
    xt::xtensor<int, 2> _corr_pair;
    float _maxlag;
    int _smoothspect_N;
    bool _flag;
    int _flag_gap;
    int _threads;

  public:
    int _npts;
    int _maxlag_npts; 
    int _fft_npts; // number of points for each channel in the fftdata
    int _channel_num; // number of channels in the fftdata 
    int _win_num; // number of windows for each channel in the fftdata
    int _pair_num; // number of pairs in the corr_pair
    int _maxlag_start; // start index of maxlag
    int _maxlag_end; // end index of maxlag
    T3 _output_data; // output data

    // constructor
    CorrClass(T2& fft_data, float dt, std::string corr_method, xt::xtensor<int, 2> corr_pair, float maxlag, int smoothspect_N, bool flag, int flag_gap, int threads);

    // copy constructor
    // CorrClass(CorrClass & A);

    // destructor
    ~CorrClass();

    // corr  deconv, coherence
    xt::xarray<std::complex<T1>>  xcorr(int pair_id);
    xt::xarray<std::complex<T1>>  deconv(int pair_id);
    xt::xarray<std::complex<T1>>  coherency(int pair_id);

    // irfft
    xt::xarray<T1> irfft(xt::xarray<std::complex<T1>>& fft_corr_data);

    // cut with maxlag
    xt::xarray<T1> cut_maxlag(xt::xarray<T1>& corr_data);

    // run workflow
    void run();

    // return
    T3& get_corrdata() {return _output_data;};
};





// **************************************************************************
// *                        constructor
// **************************************************************************
// parameter constructor                                                                                                                                                      
template <typename T1, typename T2, typename T3> 
CorrClass<T1, T2, T3>::CorrClass(T2& fft_data, float dt, std::string corr_method, xt::xtensor<int, 2> corr_pair, float maxlag, int smoothspect_N, bool flag, int flag_gap, int threads) 
{
    _fft_data = fft_data;
    _dt = dt;
    _corr_method = corr_method;
    _corr_pair = corr_pair;
    _maxlag = maxlag;
    _smoothspect_N = smoothspect_N;
    _flag = flag;
    _flag_gap = flag_gap;
    _threads = threads;

    _pair_num = corr_pair.shape(0);
    _channel_num = fft_data.shape(0);
    _win_num = fft_data.shape(1);
    _fft_npts = fft_data.shape(2);  
    _npts = 2*(_fft_npts-1);

    xt::xarray<T1> t = xt::arange(-_fft_npts+1, _fft_npts) * _dt;
    auto ind = xt::where(xt::abs(t) <= _maxlag)[0];
    _maxlag_npts = ind.size();
    _maxlag_start = ind[0]; 
    _maxlag_end = ind[ind.size()-1]+1;
}


// virtual destructor                                                                                                                                                       
template <typename T1, typename T2, typename T3> 
CorrClass<T1, T2, T3>::~CorrClass() {}






// **************************************************************************
// *                        corr
// **************************************************************************
// xcorr
template <typename T1, typename T2, typename T3> 
xt::xarray<std::complex<T1>> CorrClass<T1, T2, T3>::xcorr(int pair_id)
{
    int source_id = _corr_pair(pair_id,0);
    int receiver_id = _corr_pair(pair_id,1);
    xt::xarray<std::complex<T1>> fft_corr_data = xt::conj(xt::view(_fft_data, source_id, xt::all(), xt::all())) * xt::view(_fft_data, receiver_id, xt::all(), xt::all());

    return fft_corr_data;
}


// deconv
template <typename T1, typename T2, typename T3>
xt::xarray<std::complex<T1>> CorrClass<T1, T2, T3>::deconv(int pair_id)
{
    std::complex<T1> _i0 {0, 0};
    int source_id = _corr_pair(pair_id,0);
    int receiver_id = _corr_pair(pair_id,1);  

    auto rr_data = xt::view(_fft_data, receiver_id, xt::all(), xt::all());
    xt::xarray<std::complex<T1>> ss_data = xt::empty<std::complex<T1>>({_win_num, _fft_npts});

    for (int i=0; i<_win_num; i++){
        auto ss = xt::view(_fft_data, source_id, i, xt::all());
        auto temp = CC::moving_ave<T1>(xt::abs(ss), _smoothspect_N);
        xt::view(ss_data, i, xt::all()) = ss/temp;
    }

    // convert nan to 0+0j
    for(auto it=ss_data.begin(); it!=ss_data.end(); ++it){
        if (*it != *it){
            *it = _i0;
        }
    }

    xt::xarray<std::complex<T1>> fft_corr_data = xt::conj(ss_data) * rr_data;

  return fft_corr_data;
}


// coherency
template <typename T1, typename T2, typename T3>
xt::xarray<std::complex<T1>> CorrClass<T1, T2, T3>::coherency(int pair_id)
{
    std::complex<T1> _i0 {0, 0};
    int source_id = _corr_pair(pair_id,0);
    int receiver_id = _corr_pair(pair_id,1);   

    xt::xarray<std::complex<T1>> ss_data = xt::empty<std::complex<T1>>({_win_num, _fft_npts});
    xt::xarray<std::complex<T1>> rr_data = xt::empty<std::complex<T1>>({_win_num, _fft_npts});

    for (int i=0; i<_win_num; i++){
        auto ss = xt::view(_fft_data, source_id, i, xt::all());
        auto rr = xt::view(_fft_data, receiver_id, i, xt::all());
        auto ss_temp = CC::moving_ave<T1>(xt::abs(ss), _smoothspect_N);
        auto rr_temp = CC::moving_ave<T1>(xt::abs(rr), _smoothspect_N);
        xt::view(ss_data, i, xt::all()) = ss / ss_temp;
        xt::view(rr_data, i, xt::all()) = rr / rr_temp;
    }

    // convert nan to 0+0j
    for(auto it=ss_data.begin(); it!=ss_data.end(); ++it){
        if (*it!=*it){
            *it = _i0;
        }
    }

    // convert nan to 0+0j
    for(auto it=rr_data.begin(); it!=rr_data.end(); ++it){
        if (*it!=*it){
            *it = _i0;
        }
    }

    xt::xarray<std::complex<T1>> fft_corr_data = xt::conj(ss_data) * rr_data;

    return fft_corr_data;
}





// **************************************************************************
// *                        irfft
// **************************************************************************
template <typename T1, typename T2, typename T3>
xt::xarray<T1> CorrClass<T1, T2, T3>::irfft(xt::xarray<std::complex<T1>>& fft_corr_data)
{
    xt::xarray<T1> corr_data = xt::empty<T1>({_win_num, _npts});

    for (int i=0; i<_win_num; i++){
        auto ss = xt::view(fft_corr_data, i, xt::all());
        xt::view(corr_data, i, xt::all()) = xt::roll(xt::fftw::irfft(xt::eval(ss)), int(_fft_npts-1));
    }

    return corr_data;
}




// **************************************************************************
// *                        cut_maxlag
// **************************************************************************
template <typename T1, typename T2, typename T3>
xt::xarray<T1> CorrClass<T1, T2, T3>::cut_maxlag(xt::xarray<T1>& corr_data)
{
    auto maxlag_corr_data = xt::view(corr_data, xt::all(), xt::range(_maxlag_start, _maxlag_end));

    return maxlag_corr_data;
}





// **************************************************************************
// *                        run
// **************************************************************************
template <typename T1, typename T2, typename T3> 
void CorrClass<T1, T2, T3>::run()
{
  auto start_time = std::chrono::high_resolution_clock::now();
  auto end_time = start_time;
  _output_data = xt::empty<T1>({_pair_num, _win_num, _maxlag_npts});
  
#ifdef _OPENMP
  omp_set_num_threads(_threads);
  std::cout << "Start corr with " << _threads  << " threads using openmp..." << std::endl;
  #pragma omp parallel for
#else
  std::cout << "Start corr with " << 1 << " threads without openmp..." << std::endl;
#endif

  for(int i = 0; i < _pair_num; i++){
// flag
#ifdef _OPENMP
    if (_flag && omp_get_thread_num()==0 && i%_flag_gap==0) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << omp_get_thread_num() << " starts to process pair " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#else
    if (_flag && i%_flag_gap==0) {
      end_time = std::chrono::high_resolution_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
      std::cout << "Thread " << 0 << " starts to process pair " << i << " || " << elapsed_time.count() << "s" << std::endl;
    }
#endif

    // corr
    xt::xarray<std::complex<T1>> fft_corr_data;
    if (_corr_method == "xcorr") {
      fft_corr_data = xcorr(i);
    } 
    else if (_corr_method == "deconv"){
      fft_corr_data = deconv(i);
    } 
    else if (_corr_method == "coherency") {
      fft_corr_data = coherency(i);
    } 
    else {
      std::cout << "error: please input the correct corr method." << std::endl;
    }

    // irfft
    xt::xarray<T1> corr_data = irfft(fft_corr_data);

    // cut maxlag
    xt::xarray<T1> maxlag_corr_data = cut_maxlag(corr_data);
    
    // fill with maxlag-data
    xt::view(_output_data, i, xt::all()) = maxlag_corr_data;
  }

  end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
  std::cout << "End corr with total time " << elapsed_time.count() << "s" <<std::endl;

}



// **************************************************************************
// *                        end
// **************************************************************************

}

#endif