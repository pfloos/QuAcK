#! /bin/bash

MOL="H2"
BASIS="VDZ"

for R in $(seq 0.5 0.1 0.6)
do
  echo "# nAt nEla nElb nCore nRyd"           > examples/molecule.$MOL
  echo " 2   1     1   0     0"              >> examples/molecule.$MOL
  echo "# Znuc   x            y           z" >> examples/molecule.$MOL
  echo "  H     0.           0.         0."  >> examples/molecule.$MOL
  echo "  H     0.           0.         $R"  >> examples/molecule.$MOL
  ./GoDuck $MOL $BASIS > ${MOL}_${BASIS}_${R}.out
done

