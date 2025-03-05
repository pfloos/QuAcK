#! /bin/bash

cp ./methods.test ../input/methods
cp ./options.test ../input/options
cd ..
python3 PyDuck.py -x N2 -b sto-3g -c 0 -m 1
