#! /bin/bash

for i in *.out
do

  echo 
  echo '***********************'
  echo '***  ' $i '  ***'
  echo '***********************'
  echo 

  echo 
  echo '--- HF information ---'
  echo 

  grep "Hartree-Fock energy"       $i
  grep "HF HOMO      energy (eV):" $i
  grep "HF LUMO      energy (eV):" $i
  grep "HF HOMO-LUMO gap    (eV):" $i

  echo 
  echo '--- WFT information ---'
  echo 

  grep "Ec(MP2) =" $i
  grep "Ec(CCSD) =" $i
  grep "Ec(CCSD(T)) =" $i

  echo 
  echo '--- CIS excitation energy (singlet & triplet) ---'
  echo 
  grep "|     1 |" $i | head -1 | cut -f4 -d"|" 
  grep "|     1 |" $i | head -2 | cut -f4 -d"|" | tail -1

  echo 
  echo '--- TDHF excitation energy (singlet & triplet) ---'
  echo 
  grep "|     1 |" $i | head -3 | cut -f4 -d"|" | tail -1
  grep "|     1 |" $i | head -4 | cut -f4 -d"|" | tail -1

  echo 
  echo '--- GF2 information ---'
  echo 

  grep "GF2  HOMO      energy (eV):" $i | tail -1 
  grep "GF2  LUMO      energy (eV):" $i | tail -1 
  grep "GF2  HOMO-LUMO gap    (eV):" $i | tail -1 

  echo 
  echo '--- G0W0 information ---'
  echo 

  grep "G0W0 HOMO      energy (eV):" $i
  grep "G0W0 LUMO      energy (eV):" $i
  grep "G0W0 HOMO-LUMO gap    (eV):" $i
  grep "G0W0 RPA total energy     =" $i
  grep "G0W0 GM total energy      =" $i

  echo 
  echo '--- BSE@G0W0 excitation energy (singlet & triplet) ---'
  echo 
  grep "|     1 |" $i | head -6 | cut -f4 -d"|" | tail -1
  grep "|     1 |" $i | head -7 | cut -f4 -d"|" | tail -1

  echo 
  echo '--- evGW information ---'
  echo 

  grep "evGW HOMO      energy (eV):" $i | tail -1
  grep "evGW LUMO      energy (eV):" $i | tail -1
  grep "evGW HOMO-LUMO gap    (eV):" $i | tail -1
  grep "evGW RPA total energy     =" $i | tail -1
  grep "evGW GM total energy      =" $i | tail -1

  echo 
  echo '--- BSE@evGW excitation energy (singlet & triplet) ---'
  echo 
  grep "|     1 |" $i | tail -2 | head -1 | cut -f4 -d"|" 
  grep "|     1 |" $i | tail -1 | cut -f4 -d"|" 

  echo 
  echo '*** DONE ***'
  echo 

done

