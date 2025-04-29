#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000
ulimit -s unlimited
source ~/.bashrc
cp ./methods.test ../input/methods
cp ./options.test ../input/options
basis="aug-cc-pvtz+3s3p3d_N2"
molecule="N2X"
cp  ../cap_integrals/$basis ../int/CAP.dat
cd ..
python3.11 PyDuck.py -x $molecule -b $basis -c 0 -m 1 --use_local_basis
