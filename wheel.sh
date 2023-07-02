#!/bin/bash

# Rename the wheel file to match the currently active python version
wheel_file=$(find dist/ -type f -name '*.whl')
echo "Renaming ${wheel_file}"
python_ver=$(python -c 'import platform; print("".join(platform.python_version_tuple()[:2]))')
echo "Using Python version: ${python_ver}"
new_wheel_file=$(echo ${wheel_file} | sed "s/-cp[0-9]*/-cp${python_ver}/g")
echo "New wheel file name: ${new_wheel_file}"
mv ${wheel_file} ${new_wheel_file}