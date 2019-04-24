#! /bin/bash

Lmax=6
Mmax=6
rs=$1

if [ $# != 1 ]
then
  echo "Please, specify rs value" 
else

  echo "------------------------"
  echo "Maxmium L value = " $Lmax
  echo "Maxmium M value = " $Mmax
  echo "------------------------"
  echo
  
  for (( L=0 ; L<=$Lmax ; L++ )) ; do

    ne=$(bc -l <<< "(2*($L+1)*($L+1))")
    echo 
    echo "------------------------"
    echo "Number of electrons = " $ne
    echo "------------------------"
    echo 

    for (( M=$L+1 ; M<=$Mmax ; M++ )) ; do

      nb=$(bc -l <<< "(($M+1)*($M+1))")
      echo "Number of basis functions = " $nb
      echo -e "# rs \n" $rs > input/sph
      ./GoSph $ne $M > Sph_${ne}_${M}.out 
      grep "Total CPU time for QuAcK =" Sph_${ne}_${M}.out 

    done

  done

fi

