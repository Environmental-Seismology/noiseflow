# NoiseFlow


## Prerequisites

NoiseFlow now supports `Clang` and `GCC` compiler in MacOS and Linux system separately, and all installation processes are under the `conda` environment, and we recommend to use miniconda. Make sure to install the following pre-packages before installing noiseflow:


If you use `Clang` in Mac, please install `OpenMP` via `brew` as following:

```bash
brew install openmp
```

And use `pip` and `conda` to install the following packages:

```bash
pip install joblib

conda install -c conda-forge numpy scipy matplotlib 
conda install -c conda-forge obspy
conda install -c conda-forge fftw
conda install -c conda-forge pybind11
conda install -c conda-forge xtensor xsimd xtl xtensor-blas xtensor-python
conda install -c conda-forge xtensor-fftw  #(usually failed at most time)  
```

The `xtensor-fftw` and `KFR` need to be installed from source, first download them:


```bash
git clone https://github.com/OUCyf/noiseflow.git
cd noiseflow
git submodule init
git submodule update
```



Note the `xtensor-fftw` do not support M1 chip, and if it is failed to install via conda, you can install it from source into conda environment as `$CONDA_PREFIX`

```bash
cd ./extern/xtensor-fftw
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
make install
```



The `KFR` package is C++ DSP framework, should be installed in `./extern/kfr` from source

```bash
cd ./extern/kfr
mkdir build && cd build
cmake ..
make install
```




## Installation

Now you can install `NoiseFlow`. If you use `MacOS`, please make sure to use Clang as the complier

```bash
export CXX=clang++
# unset CXX
python setup.py install
```

If you use `Linux`, please use GCC as the compiler

```bash
export CXX=g++-13
python setup.py install
```


If you use `HPC` with `module` tool, you can use both Clang and GCC, for example using NOTS in Rice University.

```bash
# use gcc
module load GCC/13.1.0
export CXX=g++
python setup.py install
INCLUDE_CMAKE_EXTENSION=1 pip install .

# use clang
module load GCCcore/11.2.0
module load Clang/13.0.1
export CXX=clang++
python setup.py install
```

```bash
conda install -c conda-forge stockwell
```



##
Noiseflow is dual-licensed, available under both commercial and apache license.

If you want to use noiseflow in a commercial product or a closed-source project, you need to purchase a Commercial License.
