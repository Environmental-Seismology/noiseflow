# python modules
import os
import json

env = os.environ.get('CONDA_DEFAULT_ENV')
config_path = os.path.abspath(os.path.expanduser(f'~/.noiseflow_config_{env}.json'))    

try:
    from noiseflow.lib import cc_share
    from noiseflow.lib import signal_share
    NOISEFLOW_USE_CPP = True
except:
    NOISEFLOW_USE_CPP = False

compile_time_env = {"NOISEFLOW_USE_CPP": NOISEFLOW_USE_CPP}
with open(config_path, 'w') as f:
    json.dump(compile_time_env, f)


from noiseflow.cc.wrapper import rfft, corr, stack
from noiseflow.cc.rfftdata import RFFTData_Class
from noiseflow.cc.corrdata import CorrData_Class
from noiseflow.cc.stackdata import StackData_Class
from noiseflow.cc.utils_time import time_linspace, get_timestamp, get_stack_timestamp
from noiseflow.cc.utils_load import load_raw, load_rfft, load_corr, load_stack

from noiseflow.signal.rawdata import RawData_Class
from noiseflow.signal.wrapper import bandpass, bandstop, lowpass, highpass, detrend, decimate, taper

from noiseflow.client.client import downloader, downloader_https

from noiseflow.dispersion import app_dispersion, test

__all__ = ["wrapper", "dispersion", "tests"]
