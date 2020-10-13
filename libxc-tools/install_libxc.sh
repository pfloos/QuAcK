#!/bin/bash

wget https://gitlab.com/libxc/libxc/-/archive/5.0.0/libxc-5.0.0.tar.gz
tar -xzf libxc-5.0.0.tar.gz
cd libxc-5.0.0
./configure --prefix=$QUACK_ROOT
make
make check
make install
