#!/bin/bash

PKG_NAME=mpinv

echo 'Preparing required libraries, may take a while'
[ ! -d "./libs/Eigen" ] && tar -xzf libs/Eigen-3_3_8.tar.gz -C libs
[ ! -d "./libs/boost" ] && tar -xzf libs/boost-1_74_0.tar.gz -C libs

# Clean build files
rm -rf bin/*
rm -rf build/*
rm -rf env1/*
make clean


# Build wheel python library
make
python3 -m venv env1
source env1/bin/activate
pip3 install sip -q --disable-pip-version-check
sip-build

cp -r src/setup_py/ bin
cp build/${PKG_NAME}/${PKG_NAME}*.so bin/setup_py/${PKG_NAME}

python3 -m pip install --upgrade setuptools wheel
cd bin/setup_py
python3 setup.py sdist bdist_wheel

pip3 install dist/${PKG_NAME}*.tar.gz

cd ../../

python3 test/test2.py

if [ $? -eq 0 ]; then
    mkdir bin/py_module
    cp bin/setup_py/dist/* bin/py_module
    echo 'Building completed! Look into "bin/py_module/"'
    echo 'To install lib run "pip3 install --user '${PKG_NAME}'.tar.gz or .whl"'
fi