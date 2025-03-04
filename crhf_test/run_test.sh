#! /bin/bash

cp ./methods.test ../input/methods
cp ./options.test ../input/options
cd ..
python3 PyDuck.py -x N2 -b aug-cc-pVTZ -c -1 -m 2
cp input/methods.default input/methods
cp input/options.default input/options
