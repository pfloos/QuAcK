#! /bin/bash

MOL="LiF"
BASIS="cc-pvtz"
R_START=3.980
R_END=2.990
DR=0.001

for R in $(seq $R_START $DR $R_END)
do
  echo "# nAt nEla nElb nCore nRyd"           > examples/molecule.$MOL
  echo " 2   6     6   0     0"              >> examples/molecule.$MOL
  echo "# Znuc   x            y           z" >> examples/molecule.$MOL
  echo "  Li    0.           0.         0."  >> examples/molecule.$MOL
  echo "  F     0.           0.         $R"  >> examples/molecule.$MOL
  ./GoDuck $MOL $BASIS > ${MOL}_${BASIS}_${R}.out
  echo $R `./extract.sh ${MOL}_${BASIS}_${R}.out | tail -4 | head -1`
done

