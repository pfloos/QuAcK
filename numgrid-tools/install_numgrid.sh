#!/bin/bash -x

git clone https://github.com/dftlibs/numgrid.git
cd numgrid
git reset --hard 29f94b7
./setup --fc="$FC" --cc="$CC" --cxx="$CXX"
cd build
make
cp lib/libnumgrid.a $LDIR/
cp ../numgrid/numgrid.f90 $SDIR/numgrid/

