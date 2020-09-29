#! /bin/bash

MOL=$1
BASIS=$2

w_start=0.00
w_end=1.05
dw=0.05

w1=0.00

XF=$3
CF=$4

# for H
#aw1="1.49852 7.79815 25.1445"
#aw2="0.424545 -0.0382349 -0.32472"


# for He
#aw1="0.429447 0.053506 -0.339391"
#aw2="0.254939 -0.0893396 0.00765453"

# for H2
aw1="0.445525 0.0901503 -0.286898"
aw2="0.191734 -0.0364788 -0.017035"

# for Li
#aw1="0.055105 -0.00943825 -0.0267771"
#aw2="0.0359827 0.0096623 -0.0173542"

# for Li+
#aw1="0.503566, 0.137076, -0.348529"
#aw2="0.0553828, 0.00830375, -0.0234602"


# for B
#aw1="0.052676 -0.00624118 -0.000368825"
#aw2="0.0385558 -0.0015764 -0.000894297"

# for O
#aw1="-0.0187067 -0.0141017 -0.0100849"
#aw2="0.00544868 -0.0000118236 -0.000163245"

# for Al
#aw1="-0.00201219 -0.00371002 -0.00212719"
#aw2="-0.00117715 0.00188738 -0.000414761"

# for Be
#aw1="0.0663282, -0.0117682, -0.0335909"
#aw2="0.0479262, 0.00966351, -0.0208712"


DATA=${MOL}_${BASIS}_${XF}_${CF}_${w2}.dat
rm $DATA
touch $DATA

for w2 in $(seq $w_start $dw $w_end)
do
 ## w2=${w1}
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
  echo "# occupation numbers of orbitals nO and nO+1"  >> input/dft
  echo " 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo " 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo " " >> input/dft
  echo " 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo " " >> input/dft
  echo " 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo " 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " >> input/dft
  echo "# Ensemble weights: wEns(1),...,wEns(nEns-1)" >> input/dft
  echo "  ${w1} ${w2} " >> input/dft
  echo "# Ncentered ? 0 for NO " >> input/dft
  echo " 0 " >> input/dft
  echo "# Parameters for CC weight-dependent exchange functional" >> input/dft
  echo ${aw1} >> input/dft
  echo ${aw2} >> input/dft
  echo "# choice of UCC exchange coefficient : 1 for Cx1, 2 for Cx2, 3 for Cx1*Cx2" >> input/dft
  echo "2"  >> input/dft
  echo "# GOK-DFT: maxSCF thresh   DIIS n_diis guess_type ortho_type" >> input/dft
  echo "           1000    0.00001   T     5      1          1" >> input/dft
  OUTPUT=${MOL}_${BASIS}_${XF}_${CF}_${w2}.out
  ./GoXC $MOL $BASIS > ${OUTPUT}
  Ew=`grep "Ensemble energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E0=`grep "Individual energy state  1:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E1=`grep "Individual energy state  2:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  E2=`grep "Individual energy state  3:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'`
  IP=`grep "Ionization Potential"  ${OUTPUT} | grep " au" | tail -1 | cut -d":" -f 2 | sed 's/au//'`
  EA=`grep "Electronic Affinity"  ${OUTPUT} | grep " au" | tail -1 | cut -d":" -f 2 | sed 's/au//'`
  FG=`grep "Fundamental Gap"  ${OUTPUT} | grep " au" | tail -1 | cut -d":" -f 2 | sed 's/au//'`
  Ex=`grep "Exchange        energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/au//'` 
  HOMOa=`grep "HOMO a    energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/eV//'`
  LUMOa=`grep "LUMO a    energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/eV//'`
  HOMOb=`grep "HOMO a    energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/eV//'`
  LUMOb=`grep "LUMO b    energy:" ${OUTPUT} | cut -d":" -f 2 | sed 's/eV//'`

  echo $w1 $w2 $Ew $E0 $E1 $E2 $IP $EA $FG $Ex $HOMOa $LUMOa $HOMOb $LUMOb
  echo $w1 $w2 $Ew $E0 $E1 $E2 $IP $EA $FG $Ex $HOMOa $LUMOa $HOMOb $LUMOb >> ${DATA}
done

