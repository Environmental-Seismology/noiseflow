/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/


#ifndef STACK_HPP
#define STACK_HPP

#include "utils.hpp"

namespace CC {

// Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
template <typename T1, typename T2, typename T3> 
class StackClass{
  private:
    T2 _corr_data; // input data
    std::string _stack_method;
    bool _stack_all;
    float _stack_len;
    float _stack_step;
    bool _pick;
    float _median_high;
    float _median_low;
    bool _flag;
    int _flag_gap;
    int _threads;

  public:
    int _npts;
    int _pair_num; 
    int _corr_win_num;
    int _stack_win_num;
    xt::xtensor<int, 2> _win_info;
    T3 _ngood_all;
    T3 _output_data;

    // constructor
    StackClass(T2& corr_data, std::string stack_method, bool stack_all, float stack_len, float stack_step, bool pick, float median_high, float median_low, bool flag, int flag_gap, int threads);

    // copy constructor
    // StackClass(StackClass & A);

    // destructor
    ~StackClass();

    // pick --> _nindex, _ngood
    xt::xarray<int> pick(int pair_id);
    xt::xarray<T1> check_ngood(xt::xarray<int> & nindex);

    // linear, pws, robust
    xt::xarray<T1>  linear(int pair_id, xt::xarray<int> & nindex);
    xt::xarray<T1>  pws(int pair_id, xt::xarray<int> & nindex);
    xt::xarray<T1>  robust(int pair_id, xt::xarray<int> & nindex);

    // run workflow
    void run();

    // return
    T3& get_ngood() {return _ngood_all;};
    T3& get_stackdata() {return _output_data;};
    
};


// **************************************************************************
// *                        constructor
// **************************************************************************
// parameter constructor                                                                                                                                                      
template <typename T1, typename T2, typename T3> 
StackClass<T1, T2, T3>::StackClass(T2& corr_data, std::string stack_method, bool stack_all, float stack_len, float stack_step, bool pick, float median_high, float median_low, bool flag, int flag_gap, int threads) 
{ 
    _corr_data = corr_data;
    _stack_method = stack_method;
    _stack_all = stack_all;
    _stack_len = stack_len;
    _stack_step = stack_step;
    _pick = pick;
    _median_high = median_high;
    _median_low = median_low;
    _flag = flag;
    _flag_gap = flag_gap;
    _threads = threads;

    _pair_num = corr_data.shape(0);
    _corr_win_num = corr_data.shape(1);
    _npts = corr_data.shape(2);  

    // slice_window --> _stack_win_num, _win_info
    if (_stack_all){
        _win_info = xt::xtensor<int, 2> {{0, _corr_win_num}};
        _stack_win_num = 1;
    }
    else{
        _win_info = CC::slice_window(_corr_win_num, _stack_len, _stack_step);
        _stack_win_num = _win_info.shape(0);
    }
}

                                                                                                                                                       
// virtual destructor                                                                                                                                                          
template <typename T1, typename T2, typename T3> 
StackClass<T1, T2, T3>::~StackClass() {}




// **************************************************************************
// *                        pick
// **************************************************************************
template <typename T1, typename T2, typename T3> 
xt::xarray<int> StackClass<T1, T2, T3>::pick(int pair_id)
{
    xt::xarray<int> nindex;
    if (_pick){
        // SOME WARNING HERE!!!
        auto ampmax = xt::amax(xt::view(_corr_data, pair_id, xt::all(), xt::all()), 1);
        auto index = xt::where((ampmax >_median_low*xt::median(ampmax)) && (ampmax<_median_high*xt::median(ampmax)))[0]; // index: std::vector<long unsigned int>
        
        std::vector<int> index_i;
        for (auto& elem : index) {
            index_i.push_back(elem);
        }
        nindex = xt::adapt(index_i);   
    } 
    else{
        nindex = xt::arange<int>(0, _corr_win_num);
    }

    return nindex;  
}




// **************************************************************************
// *                        check_ngood
// **************************************************************************
template <typename T1, typename T2, typename T3> 
xt::xarray<T1> StackClass<T1, T2, T3>::check_ngood(xt::xarray<int> & nindex)
{
    xt::xarray<T1> ngood = xt::empty<T1>({_stack_win_num});
    if (_pick){
        for(int i=0; i<_stack_win_num; i++){
            ngood(i) = xt::sum(nindex >= _win_info(i,0) && nindex < _win_info(i,1))(0);
        }
    } 
    else {
        ngood = _stack_len*xt::ones<T1>({_stack_win_num});
    }

    return ngood;
}





// **************************************************************************
// *                        stack
// **************************************************************************
// linear
template <typename T1, typename T2, typename T3> 
xt::xarray<T1> StackClass<T1, T2, T3>::linear(int pair_id, xt::xarray<int> & nindex)
{
    xt::xarray<T1> stack_data = xt::empty<T1>({_stack_win_num, _npts});

    for(int i=0; i<_stack_win_num; i++){
        std::vector<int> pick_index;
        xt::xarray<T1> each_index = xt::arange(_win_info(i,0), _win_info(i,1));
        std::set_intersection(nindex.begin(), nindex.end(), each_index.begin(), each_index.end(), std::back_inserter(pick_index));

        auto pick_corr_data = xt::view(_corr_data, pair_id, xt::keep(pick_index), xt::all());
        if (pick_corr_data.shape(0) == 0){
            xt::view(stack_data,i, xt::all()) = xt::zeros<T1>({_npts});
        } 
        else{
            xt::view(stack_data,i, xt::all()) = xt::mean(pick_corr_data,0); // xt::sum(pick_corr_data,0) / static_cast<T1>pick_corr_data.shape(0);
        }
   }

    return stack_data;
}


// pws
template <typename T1, typename T2, typename T3> 
xt::xarray<T1> StackClass<T1, T2, T3>::pws(int pair_id, xt::xarray<int> & nindex)
{
    xt::xarray<T1> stack_data = xt::empty<T1>({_stack_win_num, _npts});

    for(int i=0; i<_stack_win_num; i++){
        std::vector<int> pick_index;
        xt::xarray<T1> each_index = xt::arange(_win_info(i,0), _win_info(i,1));
        std::set_intersection(nindex.begin(), nindex.end(), each_index.begin(), each_index.end(), std::back_inserter(pick_index));

        auto pick_corr_data = xt::view(_corr_data, pair_id, xt::keep(pick_index), xt::all());
        if (pick_corr_data.shape(0) == 0){
            xt::view(stack_data,i, xt::all()) = xt::zeros<T1>({_npts});
        } 
        else{
            xt::view(stack_data,i, xt::all()) = CC::pws<T1>(pick_corr_data);
        }
    }

    return stack_data;
}


// robust
template <typename T1, typename T2, typename T3> 
xt::xarray<T1> StackClass<T1, T2, T3>::robust(int pair_id, xt::xarray<int> & nindex)
{
    xt::xarray<T1> stack_data = xt::empty<T1>({_stack_win_num, _npts});
    
    for(int i=0; i<_stack_win_num; i++){
        std::vector<int> pick_index;
        xt::xarray<T1> each_index = xt::arange(_win_info(i,0), _win_info(i,1));
        std::set_intersection(nindex.begin(), nindex.end(), each_index.begin(), each_index.end(), std::back_inserter(pick_index));
        
        auto pick_corr_data = xt::view(_corr_data, pair_id, xt::keep(pick_index), xt::all());
        if (pick_corr_data.shape(0) == 0){
            xt::view(stack_data,i, xt::all()) = xt::zeros<T1>({_npts});
        } 
        else{
            xt::view(stack_data,i, xt::all()) = CC::robust<T1>(pick_corr_data);
        }
    }

    return stack_data;
}






// **************************************************************************
// *                        run
// **************************************************************************
template <typename T1, typename T2, typename T3> 
void StackClass<T1, T2, T3>::run()
{
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    _ngood_all = xt::empty<T1>({_pair_num, _stack_win_num});
    _output_data = xt::empty<T1>({_pair_num, _stack_win_num, _npts});
    
#ifdef _OPENMP
  omp_set_num_threads(_threads);
  std::cout << "Start stack with " << _threads  << " threads using openmp..." << std::endl;
  #pragma omp parallel for
#else
  std::cout << "Start stack with " << 1 << " threads without openmp..." << std::endl;
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

        // pick --> _nindex
        xt::xarray<int> nindex = pick(i);

        // ngood
        xt::xarray<T1> ngood = check_ngood(nindex);

        // corr
        xt::xarray<T1> stack_data;
        if (_stack_method == "linear") {
            stack_data = linear(i, nindex);
        } 
        else if (_stack_method == "pws"){
            stack_data = pws(i, nindex);
        } 
        else if (_stack_method == "robust") {
            stack_data = robust(i, nindex);
        } 
        else {
            std::cout << "error: please input the correct stack method." << std::endl;
        }
    
        // fill with maxlag-data
        xt::view(_ngood_all, i, xt::all()) = ngood;
        xt::view(_output_data, i, xt::all()) = stack_data;
    }

  end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
  std::cout << "End stack with total time " << elapsed_time.count() << "s" <<std::endl;
}




// **************************************************************************
// *                        end
// **************************************************************************

} 

#endif