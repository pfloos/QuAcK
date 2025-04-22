#!/bin/bash
#SBATCH --job-name="eta_scan"
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000

ulimit -s unlimited
source ~/.bashrc
cp ./methods.test ../input/methods
cp ./options.test ../input/options
basis="aug-cc-pvtz+3s3p3d_CO"
molecule="COX"
cp  ../cap_integrals/$basis ../int/CAP.dat
cd ..
echo "eta = "
cat ./input/eta_opt.dat
python3.11 PyDuck.py -x $molecule -b $basis -c 0 -m 1 --use_local_basis
