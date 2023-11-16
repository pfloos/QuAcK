#! /bin/bash

cp ./methods.test ../input/methods
cp ./options.test ../input/options
cd ..
python3 PyDuck.py -x He -b 6-31g -m 1
cp input/methods.default input/methods
cp input/options.default input/options
