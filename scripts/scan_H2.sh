#! /bin/bash

MOL="H2"
BASIS="cc-pvdz"
R_START=1.0
R_END=2.0
DR=0.1

for R in $(seq $R_START $DR $R_END)
do
  echo "2" > mol/${MOL}.xyz
  echo " " >> mol/${MOL}.xyz
  echo "H 0. 0. 0."  >> mol/${MOL}.xyz
  echo "H 0. 0. $(printf %f $R)"  >> mol/${MOL}.xyz
  ./GoDuck -x $MOL -b $BASIS -m 1 > ${MOL}_${BASIS}_$R.out
###  echo $R `./extract.sh ${MOL}_${BASIS}_$(printf %f $R).out | tail -4 | head -1`
done

