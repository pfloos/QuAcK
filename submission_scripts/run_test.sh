#! /bin/bash

cp ./methods.test ../input/methods
cp ./options.test ../input/options
basis=$2
molecule=$1
cp  ../cap_integrals/$basis ../int/CAP.dat
cd ..
python3.11 PyDuck.py -x $molecule -b $basis -c 0 -m 1 --use_local_basis
