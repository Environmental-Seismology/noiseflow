/***************************************************************************
 * Copyright (c) Fu Yin                                                     *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

// Note: this file is used for generated a shared library for pythonic user.


// base libraries
#include <iostream>
#include <string>

// external libraries
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#define FORCE_IMPORT_ARRAY // FORCE_IMPORT_ARRAY must be defined before including any header of xtensor-python, and must be defined only once.

// me libraries
#include "rfft.hpp"
#include "corr.hpp"
#include "stack.hpp"



// **************************************************************************
// *                        CC
// **************************************************************************
namespace CC {


// ******** float version ********
xt::pyarray<std::complex<float>> rfft_float(xt::pyarray<float>& raw_data, float dt, float cc_len, float cc_step, std::string time_norm, float clip_std, int smooth_N, std::string freq_norm, float freqmin, float freqmax, int whiten_npad, int smoothspect_N, bool flag, int flag_gap, int threads){
    // Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
    RFFTClass<float, xt::pyarray<float>, xt::pyarray<std::complex<float>>>  nfpy(raw_data, dt, cc_len, cc_step, time_norm, clip_std, smooth_N, freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N, flag, flag_gap, threads);
    nfpy.run();

    return nfpy.get_rfftdata();
}


xt::pyarray<float> corr_float(xt::pyarray<std::complex<float>>& rfft_data, float dt, std::string corr_method, xt::xtensor<int, 2> corr_pair, float maxlag, int smoothspect_N, bool flag, int flag_gap, int threads){
    // Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
    CorrClass<float, xt::pyarray<std::complex<float>>, xt::pyarray<float>>  nfpy(rfft_data, dt, corr_method, corr_pair, maxlag, smoothspect_N, flag, flag_gap, threads);
    nfpy.run();

    return nfpy.get_corrdata();
}


// ******** double version ********
xt::pyarray<std::complex<double>> rfft_double(xt::pyarray<double>& raw_data, float dt, float cc_len, float cc_step, std::string time_norm, float clip_std, int smooth_N, std::string freq_norm, float freqmin, float freqmax, int whiten_npad, int smoothspect_N, bool flag, int flag_gap, int threads){
    // Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
    RFFTClass<double, xt::pyarray<double>, xt::pyarray<std::complex<double>>>  nfpy(raw_data, dt, cc_len, cc_step, time_norm, clip_std, smooth_N, freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N, flag, flag_gap, threads);
    nfpy.run();

    return nfpy.get_rfftdata();
}


xt::pyarray<double> corr_double(xt::pyarray<std::complex<double>>& rfft_data, float dt, std::string corr_method, xt::xtensor<int, 2> corr_pair, float maxlag, int smoothspect_N, bool flag, int flag_gap, int threads){
    // Template (T1, T2, T3) = (data-type, input-array-type, output-array-type)
    CorrClass<double, xt::pyarray<std::complex<double>>, xt::pyarray<double>>  nfpy(rfft_data, dt, corr_method, corr_pair, maxlag, smoothspect_N, flag, flag_gap, threads);
    nfpy.run();

    return nfpy.get_corrdata();
}

}




// **************************************************************************
// *                        PYBIND11_MODULE
// **************************************************************************
namespace py = pybind11;

PYBIND11_MODULE(cc_share,m) 
{
    xt::import_numpy();

    m.doc() = "Cross-correlation module for NoiseFlow";
    
    m.def("rfft_float", &CC::rfft_float, R"pbdoc(Do rfft in c++)pbdoc");
    m.def("rfft_double", &CC::rfft_double, R"pbdoc(Do rfft in c++)pbdoc");

    m.def("corr_float", &CC::corr_float, R"pbdoc(Do correlation in c++)pbdoc");
    m.def("corr_double", &CC::corr_double, R"pbdoc(Do correlation in c++)pbdoc");

    py::class_<CC::StackClass<float, xt::pyarray<float>, xt::pyarray<float>>>(m, "StackClass_float")
        .def(py::init<xt::pyarray<float>&, std::string, bool, float, float, bool, float, float, bool, int, int>())
        .def("run", &CC::StackClass<float, xt::pyarray<float>, xt::pyarray<float>>::run)
        .def("get_ngood", &CC::StackClass<float, xt::pyarray<float>, xt::pyarray<float>>::get_ngood)
        .def("get_stackdata", &CC::StackClass<float, xt::pyarray<float>, xt::pyarray<float>>::get_stackdata);

    py::class_<CC::StackClass<double, xt::pyarray<double>, xt::pyarray<double>>>(m, "StackClass_double")
        .def(py::init<xt::pyarray<double>&, std::string, bool, float, float, bool, float, float, bool, int, int>())
        .def("run", &CC::StackClass<double, xt::pyarray<double>, xt::pyarray<double>>::run)
        .def("get_ngood", &CC::StackClass<double, xt::pyarray<double>, xt::pyarray<double>>::get_ngood)
        .def("get_stackdata", &CC::StackClass<double, xt::pyarray<double>, xt::pyarray<double>>::get_stackdata);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

}
