###################################################################  
#   Optimization Parameters
################################################################### 
cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -g -fPIC -O3")
set(XTENSOR_USE_XSIMD 1) # set xtensor to use xsimd for SIMD optimization
set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX}) # set conda env path, so that cmake can find those libraries in conda env




###################################################################  
#   Project Info & Compiler Detection
###################################################################  
project(noiseflow_share VERSION 0.1.0 LANGUAGES CXX)
message("Check C++ compiler: ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (MINGW)
        message("Using MinGW-w64")
        set(COMPILER "MinGW-w64")
    else()
        message("Using GCC")
        set(COMPILER "GCC")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message("Using Clang")
    set(COMPILER "Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message("Using AppleClang")
    set(COMPILER "AppleClang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("Using MSVC")
    set(COMPILER "MSVC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("Using Intel Compiler")
    set(COMPILER "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
    message("Using Intel oneAPI Compiler")
    set(COMPILER "Intel oneAPI")
else()
    message(FATAL_ERROR "Unrecognized compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()




###################################################################  
#   Find Libraries
###################################################################  
### [a]. use find_package for numpy
find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
### [b]. use execute_process for numpy
# execute_process(
#   COMMAND "${PYTHON_EXECUTABLE}"
#   -c "import numpy; print(numpy.get_include())"
#   OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )


### pybind11 xtensor xtl xsimd xtensor-blas xtensor-python
find_package(pybind11 REQUIRED CONFIG)
find_package(xtensor REQUIRED)
find_package(xtl REQUIRED)
find_package(xsimd REQUIRED)
find_package(xtensor-blas REQUIRED)
find_package(xtensor-python REQUIRED)


### xtensor-fftw
set(FFTW_USE_FLOAT ON CACHE BOOL "Use float precision") # set fftw3 to use float precision
set(FFTW_USE_DOUBLE ON CACHE BOOL "Use double precision") # set fftw3 to use double precision
set(FFTW_USE_LONG_DOUBLE OFF CACHE BOOL "Use long double precision") # set fftw3 to use long double precision
# find_package(xtensor-fftw REQUIRED) # this command must be after set(FFTW_USE_FLOAT ON CACHE BOOL "Use double precision")
include_directories("./extern/xtensor-fftw/include")


### FFTW3 using FindFFTW.cmake (note this command will seach some default path to find fftw3, for example the conda env path)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/extern/xtensor-fftw/cmake/Modules/findFFTW") # use FindFFTW.cmake to find fftw3
find_package(FFTW REQUIRED COMPONENTS ${REQUIRE_FLOAT_LIB} ${REQUIRE_DOUBLE_LIB} ${REQUIRE_LONG_DOUBLE_LIB})
if(FFTW_FOUND)
    message(STATUS "FFTW Found!")
endif()


### OpenBLAS and BLAS/LAPACK, which is xtensor-blas' dependencies
add_definitions(-DHAVE_CBLAS=1)
# if (WIN32)
#     find_package(OpenBLAS REQUIRED)
#     set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
# else()
#     find_package(BLAS REQUIRED)
#     find_package(LAPACK REQUIRED)
# endif()
# message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
# message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


### OpenMP (clang doesn't bring openmp, you need to install it by brew, and find_package can't find it, you need to set the PATH as below. refer to https://www.kancloud.cn/csyangbinbin/cmake-cookbook1/2157938  https://iscinumpy.gitlab.io/post/omp-on-high-sierra/)
if(APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
    message(STATUS "Using AppleClang in mac will set brew lib path for OpenMP")
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE BREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${BREW_PREFIX}/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${BREW_PREFIX}/opt/libomp/lib")
    find_package(OpenMP REQUIRED)
else()
    find_package(OpenMP REQUIRED)
endif()


### KFR 
include_directories("./extern/kfr/include")




###################################################################  
#   Target -- Python Package Module
################################################################### 
pybind11_add_module(cc_share ${CMAKE_SOURCE_DIR}/noiseflow/cc/src/pybind11.cpp)
target_include_directories(cc_share PRIVATE ${CMAKE_SOURCE_DIR}/noiseflow/cc/src ${CMAKE_SOURCE_DIR}/noiseflow/cc/include)
if(OPENMP_FOUND)
    target_link_libraries(cc_share PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} pybind11::module xtensor xtensor::optimize xtensor::use_xsimd xtensor-python Python::NumPy OpenMP::OpenMP_CXX) 
else()
    target_link_libraries(cc_share PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} pybind11::module xtensor xtensor::optimize xtensor::use_xsimd xtensor-python)
endif()


pybind11_add_module(signal_share ${CMAKE_SOURCE_DIR}/noiseflow/signal/src/pybind11.cpp)
target_include_directories(signal_share PRIVATE ${CMAKE_SOURCE_DIR}/noiseflow/signal/src ${CMAKE_SOURCE_DIR}/noiseflow/signal/include)
if(OPENMP_FOUND)
    target_link_libraries(signal_share PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} pybind11::module xtensor xtensor::optimize xtensor::use_xsimd xtensor-python Python::NumPy OpenMP::OpenMP_CXX)
else()
    target_link_libraries(signal_share PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} pybind11::module xtensor xtensor::optimize xtensor::use_xsimd xtensor-python)
endif()


### setup.py will install the module to the python site-packages
# install(TARGETS cc_share DESTINATION lib)
# install(TARGETS signal_share DESTINATION lib)
