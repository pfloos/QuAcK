#! /bin/bash

cp ./methods.test ../input/methods
cp ./options.test ../input/options
cd ..
python3 PyDuck.py -x He -b 6-31g -m 1
