#!/bin/bash



${PYTHON} -m pip install poetry-core@https://github.com/python-poetry/poetry-core/archive/refs/tags/1.6.1.zip



git clone https://github.com/xtensor-stack/xtensor-fftw extern/xtensor-fftw
git clone https://github.com/kfrlib/kfr extern/kfr

# cd ./extern/xtensor-fftw 
# mkdir build && cd build 
# cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DFFTW_USE_LONG_DOUBLE=OFF 
# make install
# cd ../../..

# cd ./extern/kfr 
# mkdir build && cd build 
# cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX 
# make install
# cd ../../..



NOISEFLOW_USE_CPP=1 ${PYTHON} -m pip install . --no-deps --ignore-installed -vvv

