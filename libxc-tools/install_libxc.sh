#!/bin/bash

###git clone --single-branch --branch master https://gitlab.com/libxc/libxc.git
tar -xzvf libxc-5.0.0.tar.gz
cd libxc-5.0.0
###autoreconf -i
./configure --prefix=$QUACK_ROOT
make
make check
make install
