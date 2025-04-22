#!/bin/bash
#SBATCH --job-name="eta_scan"
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32000

ulimit -s unlimited
source ~/.bashrc
echo "$eta" > "$QUACK_ROOT/input/eta_opt.dat"
workdir=workdir$eta
mkdir $workdir
cp -r input $workdir/input
cp -r basis $workdir/basis
cp -r cap_integrals $workdir/cap_integrals
cap -r mol $workdir/mol
cap -r cap_data $workdir/cap_data
basis="aug-cc-pvtz+3s3p3d_CO"
molecule="COX"
cp  ../cap_integrals/$basis ../int/CAP.dat
cd ..
echo "eta = "
cat ./input/eta_opt.dat
python3.11 PyDuck.py -x $molecule -b $basis -c 0 -m 1 --use_local_basis --working_dir $working_dir
