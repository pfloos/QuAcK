#! /bin/bash

MOL="N2"
BASIS="cc-pvqz"
R_START=1.7
R_END=2.6
DR=0.1

for R in $(seq $R_START $DR $R_END)
do
  echo "# nAt nEla nElb nCore nRyd"           > examples/molecule.$MOL
  echo " 2   7     7   4     0"              >> examples/molecule.$MOL
  echo "# Znuc   x            y           z" >> examples/molecule.$MOL
  echo "  N     0.           0.         0."  >> examples/molecule.$MOL
  echo "  N     0.           0.         $R"  >> examples/molecule.$MOL
  ./GoDuck $MOL $BASIS > ${MOL}_${BASIS}_FC_${R}.out
  echo $R `./extract.sh ${MOL}_${BASIS}_FC_${R}.out | tail -4 | head -1`
done

