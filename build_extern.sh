#!/bin/bash

git submodule update --init --recursive

cd ./extern/xtensor-fftw 
mkdir build && cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DFFTW_USE_LONG_DOUBLE=OFF 
make install
cd ../../..

cd ./extern/kfr 
mkdir build && cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX 
make install
cd ../../..
