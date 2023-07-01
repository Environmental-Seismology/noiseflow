#!/bin/bash

# git submodule update --init --recursive

# cd ./external/xtensor-fftw
# cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFI}
# make install
# cd ../..

# cd ./external/kfr
# mkdir build && cd build
# cmake .. -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}
# make install
# cd ../../..


# # brew install omplib
# conda install -y -c conda-forge obspy

# conda install -y -c conda-forge fftw
# conda install -y -c conda-forge pybind11
# conda install -y -c conda-forge xtensor xsimd xtl xtensor-blas xtensor-python
${PYTHON} -m pip install obspy
${PYTHON} -m pip install poetry-core@https://github.com/python-poetry/poetry-core/archive/refs/tags/1.3.2.zip

# poetry install --only main

${PYTHON} -m pip install . -vvv
# NOISEFLOW_USE_CPP=1 
# poetry build
# ${PYTHON} -m pip install dist/*.whl
