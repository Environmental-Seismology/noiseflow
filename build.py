import os
import re
import sys
import json
import time
import platform
import subprocess

from setuptools import Extension
from distutils.version import LooseVersion
from setuptools.command.build_ext import build_ext



class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)


def check_env():
    CXX = os.getenv('CXX', '')
    NOISEFLOW_USE_CPP = (os.getenv('NOISEFLOW_USE_CPP', '') == '1')
    compile_time_env = {"NOISEFLOW_USE_CPP": NOISEFLOW_USE_CPP, "CXX": CXX}
    
    config_path = os.path.abspath(os.path.expanduser('~/.noiseflow_config.json'))

    if os.path.exists(config_path):
        os.remove(config_path)
    
    with open(config_path, 'w') as f:
        json.dump(compile_time_env, f)

    ext_modules = []
    if NOISEFLOW_USE_CPP:
        ext_modules.append(CMakeExtension('noiseflow/lib/'))
        try:
            import numpy
        except:
            raise Exception(
                "please install mpi4py first to enable the MPI mode!")

    return ext_modules



def build(setup_kwargs):
    """
    This is a callback for poetry used to hook in our extensions.
    """
    setup_kwargs.update(
        {
            "ext_modules": check_env(),
            "cmdclass": {"build_ext": CMakeBuild},
        }
    )


# poetry build &&  pip install ./dist/noiseflow-0.0.1.tar.gz
