import re
import os
import argparse
import platform
import subprocess
"""
Usage:
    python build_extern.py -c [source|conda]
"""

### [build-system]
cmake = '^3.18.0'
fftw = '^3.3.10'
pybind11 = '^2.10.4'
xtensor = '^0.24.6'
xsimd = '8.0.5' # when the version is larger than 8.0.5, it will got error in compiling stage.
xtl = '^0.7.5'
xtensor_blas = '^0.20.0'
xtensor_python = '^0.26.1'
xtensor_fftw = '^0.2.6' # not show in conda env
kfr = '^5.0.2' # not show in conda env

### [optional-for-mac]
libomp = '^16.0.6' # via brew




#############################################
### functions
#############################################
def get_conda_package_info(package_name):
    try:
        output = subprocess.check_output(['conda', 'list', package_name])
        output = output.decode('utf-8')
        lines = output.split('\n')
        for line in lines:
            if line.startswith(package_name):
                package_info = re.split(r'\s+', line.strip())
                try:
                    version = package_info[1]
                except:
                    version = None
                try:
                    channel = package_info[3]
                except:
                    channel = None
                return version, channel
    except subprocess.CalledProcessError:
        pass

    return None, None


def get_channel():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', 
                    '--channel', 
                    help='specify the download source', 
                    required=False)

    args = parser.parse_args()
    if args.channel == None:
        channel = 'conda'
    else:
        channel = args.channel

    return channel


def check_libomp_installed():
    try:
        output = subprocess.check_output("brew list", shell=True)
        if 'libomp' in output.decode('utf-8'):
            return True
        else:
            return False
    except subprocess.CalledProcessError:
        print("An error occurred while checking for libomp")
        return False
    



#############################################
#### 1. get -c [channel]
#############################################
channel = get_channel()




#############################################
#### 2. install libomp via brew only for mac
#############################################
if platform.system() == "Darwin":
    if check_libomp_installed() == False:
        try:
            os.system('brew install libomp')
        except:
            print("An error occurred while installing libomp via brew")
    else:
        print("libomp has already been installed via brew")




#############################################
#### 3. must be from conda-forge
#############################################
cmake_version, cmake_channel = get_conda_package_info('cmake')
print(f'cmake version is {cmake_version} from {cmake_channel}')
if cmake_version == None:
    os.system('conda install -y -c conda-forge cmake'  )
    
fftw_version, fftw_channel = get_conda_package_info('fftw')
print(f'fftw version is {fftw_version} from {fftw_channel}')
if fftw_version == None:
    os.system('conda install -y -c conda-forge fftw'  )

pybind11_version, pybind11_channel = get_conda_package_info('pybind11')
print(f'pybind11 version is {pybind11_version} from {pybind11_channel}')
if pybind11_version == None:
    os.system('conda install -y -c conda-forge pybind11'  )




#############################################
#### 4. conda-forge or source
#############################################
if channel == 'conda':
    xtensor_version, xtensor_channel = get_conda_package_info('xtensor')
    print(f'xtensor version is {xtensor_version} from {xtensor_channel}')
    if xtensor_version == None:
        os.system('conda install -y -c conda-forge xtensor'  )
    
    xsimd_version, xsimd_channel = get_conda_package_info('xsimd')
    print(f'xsimd version is {xsimd_version} from {xsimd_channel}')
    if xsimd_version == None:
        os.system('conda install -y -c conda-forge xsimd==8.0.5'  )
    
    xtl_version, xtl_channel = get_conda_package_info('xtl')
    print(f'xtl version is {xtl_version} from {xtl_channel}')
    if xtl_version == None:
        os.system('conda install -y -c conda-forge xtl' )
    
    xtensor_blas_version, xtensor_blas_channel = get_conda_package_info('xtensor-blas')
    print(f'xtensor-blas version is {xtensor_blas_version} from {xtensor_blas_channel}')
    if xtensor_blas_version == None:
        os.system('conda install -y -c conda-forge xtensor-blas' )
    
    xtensor_python_version, xtensor_python_channel = get_conda_package_info('xtensor-python')
    print(f'xtensor-python version is {xtensor_python_version} from {xtensor_python_channel}')
    if xtensor_python_version == None:
        os.system('conda install -y -c conda-forge xtensor-python' )

elif channel == 'source':
    os.system('git submodule update --init --recursive' )

    os.system("git submodule foreach --recursive 'git clean -xfd'" )
    
    os.system("cd ./extern/xtensor \
                   && mkdir build && cd build \
                   && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                   && make install" )

    os.system("cd ./extern/xtl \
                   && mkdir build && cd build \
                   && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                   && make install" )
    
    os.system("cd ./extern/xsimd \
                   && mkdir build && cd build \
                   && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                   && make install" )

    os.system("cd ./extern/xtensor-blas \
                   && mkdir build && cd build \
                   && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                   && make install", 
                   shell=True, 
                   check=True)
    
    os.system("cd ./extern/xtensor-python \
                   && mkdir build && cd build \
                   && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                   && make install" )

else:
    raise ValueError('The parameter `-c` [channel] must be `source` or `conda`')




#############################################
#### 5. must be from source
#############################################
os.system('git submodule update --init --remote extern/kfr')

os.system('git submodule update --init --remote extern/xtensor-fftw' )

os.system("git submodule foreach --recursive 'git clean -xfd'" )

os.system("cd ./extern/xtensor-fftw \
                && mkdir build && cd build \
                && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DFFTW_USE_LONG_DOUBLE=OFF \
                && make install" )

os.system("cd ./extern/kfr \
                && mkdir build && cd build \
                && cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
                && make install" )
