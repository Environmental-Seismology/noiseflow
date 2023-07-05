#!/bin/bash



${PYTHON} -m pip install poetry-core@https://github.com/python-poetry/poetry-core/archive/refs/tags/1.6.1.zip

# ${PYTHON} -m pip install obspy

# conda env update -f environment.yml
${PYTHON} build_extern.py -c github-actions
# conda install -y -c conda-forge obspy
# ${PYTHON} -m pip install . -vvv

poetry install --only main
${PYTHON} -m pip install .
# NOISEFLOW_USE_CPP=1 poetry build
