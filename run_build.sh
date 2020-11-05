#!/bin/bash

PKG_NAME=py_mpinv2

#Prepare required libraries
[ ! -d "./libs/Eigen" ] && tar -xzf libs/Eigen-3_3_8.tar.gz -C libs
[ ! -d "./libs/boost" ] && tar -xzf libs/boost-1_74_0.tar.gz -C libs


# Clean build files
rm -rf build/*
rm -rf env1/*
make clean


# Build wheel python library
make
python3 -m venv env1
source env1/bin/activate
pip3 install sip -q --disable-pip-version-check
sip-wheel --name ${PKG_NAME}

mkdir bin/py_module
cp src/mpinv.py bin/py_module

mv ${PKG_NAME}*.whl bin/py_module

# Install lib and run python test (in env)
pip3 install mpmath -q --disable-pip-version-check
pip3 install bin/py_module/${PKG_NAME}*.whl --disable-pip-version-check

if [ $? -eq 0 ]; then
    python3 ./test/test1.py
fi

if [ $? -eq 0 ]; then
    echo 'Building completed! Look at "bin/py_module/"'
    echo 'To install lib run "pip3 install --user '${PKG_NAME}'-blabla.whl" and use "mpinv.py" as module for import'
fi