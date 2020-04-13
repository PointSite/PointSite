#!/bin/bash

#-> 1. create an conda environment 
conda remove --name pointsite_inference --all -y
conda create --name pointsite_inference python=3.6 -y
source activate pointsite_inference

#-> 2. install related library
pip install torch==1.1.0
conda install scikit-learn -y
conda install tqdm -y

#-> 3. run setup
rm -rf build/ dist/ sparseconvnet.egg-info sparseconvnet/SCN*.so
python setup.py develop
python util/setup_test.py
rm -rf build/ dist/ sparseconvnet.egg-info 

#-> 4. deactivate
source deactivate 2> /dev/null
