#! /bin/bash

MOL="H2"
BASIS="AVTZ"
R_START=0.5
R_END=3.0
DR=0.01

for R in $(seq $R_START $DR $R_END)
do
  echo "# nAt nEla nElb nCore nRyd"           > examples/molecule.$MOL
  echo " 2   1     1   0     0"              >> examples/molecule.$MOL
  echo "# Znuc   x            y           z" >> examples/molecule.$MOL
  echo "  H     0.           0.         0."  >> examples/molecule.$MOL
  echo "  H     0.           0.         $R"  >> examples/molecule.$MOL
  ./GoDuck $MOL $BASIS > ${MOL}_${BASIS}_${R}.out
  echo $R `./extract.sh ${MOL}_${BASIS}_${R}.out | tail -2 | head -1`
done

