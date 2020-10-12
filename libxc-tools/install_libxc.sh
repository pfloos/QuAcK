#!/bin/bash

git clone --single-branch --branch master https://gitlab.com/libxc/libxc.git
cd libxc
autoreconf -i
./configure --prefix=$QUACK_ROOT
make
make check
make install
