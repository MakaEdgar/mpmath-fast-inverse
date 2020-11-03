#!/bin/bash

WHL_NAME=py_mpinv2-0.1-cp36-cp36m-manylinux1_x86_64.whl


#Prepare required libraries
[ ! -d "./libs/Eigen" ] && tar -xzf libs/Eigen-3_3_8.tar.gz -C libs
[ ! -d "./libs/boost" ] && tar -xzf libs/boost-1_74_0.tar.gz -C libs


# Clean build files
rm -rf build/*
rm -rf env1/*
make clean


# Build wheel python library
make
sip-wheel --name py-mpinv2
mv ${WHL_NAME} bin


# Install lib and run python (in env)
python3 -m venv env1
source env1/bin/activate
pip3 install bin/${WHL_NAME}
python3 -ic 'from mpinv2 import mpmat;a=mpmat("./test/8x8.mpmat");print(a.version_)'
