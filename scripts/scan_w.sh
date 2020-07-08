#! /bin/bash

MOL=$1
BASIS=$2

w_start=0.0
w_end=1.0
dw=0.1

w2=0.0

XF=$3
CF=$4

aw1="0.000000 0.0000000 0.000000"
aw2="0.000000 0.0000000 0.0000000"

DATA=${MOL}_${BASIS}_${XF}_${CF}.dat
rm $DATA
touch $DATA

for w1 in $(seq $w_start $dw $w_end)
do
###  w2=${w1}
  echo "# Restricted or unrestricted KS calculation" > input/dft
  echo "  eDFT-UKS" >> input/dft
  echo "# exchange rung:" >> input/dft
  echo "# Hartree      = 0" >> input/dft
  echo "# LDA          = 1: RS51,RMFL20" >> input/dft
  echo "# GGA          = 2: RB88" >> input/dft
  echo "# Hybrid       = 4" >> input/dft
  echo "# Hartree-Fock = 666" >> input/dft
  echo "  1 $XF  " >> input/dft
  echo "# correlation rung: " >> input/dft
  echo "# Hartree      = 0" >> input/dft
  echo "# LDA          = 1: RVWN5,RMFL20" >> input/dft
  echo "# GGA          = 2: " >> input/dft
  echo "# Hybrid       = 4: " >> input/dft
  echo "# Hartree-Fock = 666" >> input/dft
  echo "  0 $CF " >> input/dft
  echo "# quadrature grid SG-n" >> input/dft
  echo "  1" >> input/dft
  echo "# Number of states in ensemble (nEns)" >> input/dft
  echo "  3" >> input/dft
  echo "# Ensemble weights: wEns(1),...,wEns(nEns-1)" >> input/dft
  echo "  ${w1} ${w2} " >> input/dft
  echo "# Parameters for CC weight-dependent exchange functional" >> input/dft
  echo ${aw1} >> input/dft
  echo ${aw2} >> input/dft
  echo "# GOK-DFT: maxSCF thresh   DIIS n_diis guess_type ortho_type" >> input/dft
  echo "           32    0.00001   T     5      1          1" >> input/dft
  OUTPUT=${MOL}_${BASIS}_${XF}_${CF}_${w1}.out
  ./GoXC $MOL $BASIS > ${OUTPUT}
  Ew=`grep "Ensemble energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E0=`grep "Individual energy state  1:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E1=`grep "Individual energy state  2:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E2=`grep "Individual energy state  3:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  IP=`grep "Ionization Potential"  ${OUTPUT} | grep " au" | tail -1 | cut -d":" -f 2 | sed 's/au//'`
  EA=`grep "Electronic Affinity"  ${OUTPUT} | grep " au" | tail -1 | cut -d":" -f 2 | sed 's/au//'`
  echo $w1 $w2 $Ew $E0 $E1 $E2 $IP $EA
  echo $w1 $w2 $Ew $E0 $E1 $E2 $IP $EA >> ${DATA}
done

