#! /bin/bash

MOL="F2"
BASIS="cc-pvtz"
R_START=2.0
R_END=3.4
DR=0.1

for R in $(seq $R_START $DR $R_END)
do
  echo "# nAt nEla nElb nCore nRyd"           > examples/molecule.$MOL
  echo " 2   9     9   0     0"              >> examples/molecule.$MOL
  echo "# Znuc   x            y           z" >> examples/molecule.$MOL
  echo "  F     0.           0.         0."  >> examples/molecule.$MOL
  echo "  F     0.           0.         $R"  >> examples/molecule.$MOL
  ./GoDuck $MOL $BASIS > ${MOL}_${BASIS}_${R}.out
  echo $R `./extract.sh ${MOL}_${BASIS}_${R}.out | tail -4 | head -1`
done

