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
#include "filter.hpp"
#include "decimate.hpp"
#include "detrend.hpp"
#include "taper.hpp"



// **************************************************************************
// *                        SIGNALL
// **************************************************************************
namespace SIGNAL {


// **************************************************************************
// *                        filter
// **************************************************************************
// bandpass
void bandpass_float(xt::pyarray<float> &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    bandpass<float, xt::pyarray<float>>(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads);
}

void bandpass_double(xt::pyarray<double> &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    bandpass<double, xt::pyarray<double>>(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads);
}


// bandstop
void bandstop_float(xt::pyarray<float> &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    bandstop<float, xt::pyarray<float>>(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads);
}

void bandstop_double(xt::pyarray<double> &data, double freqmin, double freqmax, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    bandstop<double, xt::pyarray<double>>(data, freqmin, freqmax, df, corners, zerophase, flag, flag_gap, threads);
}



// lowpass
void lowpass_float(xt::pyarray<float> &data, double freqlow, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    lowpass<float, xt::pyarray<float>>(data, freqlow, df, corners, zerophase, flag, flag_gap, threads);
}

void lowpass_double(xt::pyarray<double> &data, double freqlow, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    lowpass<double, xt::pyarray<double>>(data, freqlow, df, corners, zerophase, flag, flag_gap, threads);
}


// highpass
void highpass_float(xt::pyarray<float> &data, double freqhigh, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    highpass<float, xt::pyarray<float>>(data, freqhigh, df, corners, zerophase, flag, flag_gap, threads);
}

void highpass_double(xt::pyarray<double> &data, double freqhigh, double df, int corners, bool zerophase, bool flag, int flag_gap, int threads)
{
    highpass<double, xt::pyarray<double>>(data, freqhigh, df, corners, zerophase, flag, flag_gap, threads);
}



// **************************************************************************
// *                        DECIMATE
// **************************************************************************
xt::pyarray<float> decimate_float(xt::pyarray<float>  &data, double df, double df_new, bool flag, int flag_gap, int threads)
{
    xt::pyarray<float> results = decimate<float, xt::pyarray<float>>(data, df, df_new, flag, flag_gap, threads);

    return results;
}


xt::pyarray<double> decimate_double(xt::pyarray<double>  &data, double df, double df_new, bool flag, int flag_gap, int threads)
{
    xt::pyarray<double> results = decimate<double, xt::pyarray<double>>(data, df, df_new, flag, flag_gap, threads);

    return results;
}


// **************************************************************************
// *                        DETREND
// **************************************************************************
void detrend_float(xt::pyarray<float>  &data, std::string type, bool flag, int flag_gap, int threads)
{
    detrend<float, xt::pyarray<float>>(data, type, flag, flag_gap, threads);
}    


void detrend_double(xt::pyarray<double>  &data, std::string type, bool flag, int flag_gap, int threads)
{
    detrend<double, xt::pyarray<double>>(data, type, flag, flag_gap, threads);
}


// **************************************************************************
// *                        TAPER
// **************************************************************************
void taper_float(xt::pyarray<float> &data, double max_percentage, std::string type, std::string side, bool flag, int flag_gap, int threads) 
{
    taper<float, xt::pyarray<float>>(data, max_percentage, type, side, flag, flag_gap, threads);
}


void taper_double(xt::pyarray<double> &data, double max_percentage, std::string type, std::string side, bool flag, int flag_gap, int threads) 
{
    taper<double, xt::pyarray<double>>(data, max_percentage, type, side, flag, flag_gap, threads);
}




} // namespace SIGNAL
// **************************************************************************
// *                        PYBIND11_MODULE
// **************************************************************************
namespace py = pybind11;

PYBIND11_MODULE(signal_share,m) 
{
    xt::import_numpy();

    m.doc() = "Signal module for DasFlow";
    
    m.def("bandpass_float", &SIGNAL::bandpass_float, R"pbdoc(Do bandpass in c++)pbdoc");
    m.def("bandstop_float", &SIGNAL::bandstop_float, R"pbdoc(Do bandstop in c++)pbdoc");
    m.def("lowpass_float", &SIGNAL::lowpass_float, R"pbdoc(Do lowpass in c++)pbdoc");
    m.def("highpass_float", &SIGNAL::highpass_float, R"pbdoc(Do highpass in c++)pbdoc");
    m.def("decimate_float", &SIGNAL::decimate_float, R"pbdoc(Do decimate in c++)pbdoc");
    m.def("detrend_float", &SIGNAL::detrend_float, R"pbdoc(Do detrend in c++)pbdoc");
    m.def("taper_float", &SIGNAL::taper_float, R"pbdoc(Do taper in c++)pbdoc");

    m.def("bandpass_double", &SIGNAL::bandpass_double, R"pbdoc(Do bandpass in c++)pbdoc");
    m.def("bandstop_double", &SIGNAL::bandstop_double, R"pbdoc(Do bandstop in c++)pbdoc");
    m.def("lowpass_double", &SIGNAL::lowpass_double, R"pbdoc(Do lowpass in c++)pbdoc");
    m.def("highpass_double", &SIGNAL::highpass_double, R"pbdoc(Do highpass in c++)pbdoc");
    m.def("decimate_double", &SIGNAL::decimate_double, R"pbdoc(Do decimate in c++)pbdoc");
    m.def("detrend_double", &SIGNAL::detrend_double, R"pbdoc(Do detrend in c++)pbdoc");
    m.def("taper_double", &SIGNAL::taper_double, R"pbdoc(Do taper in c++)pbdoc");



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

}
