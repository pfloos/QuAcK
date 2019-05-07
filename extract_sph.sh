#! /bin/bash

INPUT=$1

  echo 
  echo '******************************************'
  echo '*** Extracting information of' $INPUT ' ***'
  echo '******************************************'
  echo 

  echo 
  echo '*** WFT information ***'
  grep "MP2 correlation energy" $INPUT
  grep "Ec(MP2) =" $INPUT
  grep "Ec(CCSD) =" $INPUT
  grep "Ec(CCSD(T)) =" $INPUT

  echo 
  echo '*** Gap information: HF, G0F2, GF2, G0W0 & evGW ***'
  HF=`grep "HF HOMO-LUMO gap    (eV):" $INPUT | cut -f2 -d":"`
  G0F2=`grep "GF2  HOMO-LUMO gap    (eV):" $INPUT | head -1 | cut -f2 -d":"`
  GF2=`grep "GF2  HOMO-LUMO gap    (eV):" $INPUT | tail -1 | cut -f2 -d":"`
  G0W0=`grep "G0W0 HOMO-LUMO gap    (eV):" $INPUT | cut -f2 -d":"`
  evGW=`grep "evGW HOMO-LUMO gap    (eV):" $INPUT | tail -1 | cut -f2 -d":"`

  echo -e "\t" $HF "\t" $G0F2 "\t" $GF2 "\t" $G0W0 "\t" $evGW

  echo 
  echo '*** Ec@G0W0 information: RPA, GM, BSE1 & BSE3 ***'
  RPA_G0W0=`grep "RPA@G0W0 correlation energy =" $INPUT| cut -f2 -d"="`
  GM_G0W0=`grep "GM@G0W0  correlation energy =" $INPUT| cut -f2 -d"="`
  BSE1_G0W0=`grep "BSE@G0W0 correlation energy (singlet)" $INPUT| cut -f2 -d"="`
  BSE3_G0W0=`grep "BSE@G0W0 correlation energy (triplet)" $INPUT| cut -f2 -d"="`

  echo -e "\t" $RPA_G0W0 "\t" $GM_G0W0 "\t" $BSE1_G0W0 "\t" $BSE3_G0W0

  echo 
  echo '*** Ec@evGW information: RPA, GM, BSE1 & BSE3 ***'
  RPA_evGW=`grep "RPA@evGW correlation energy =" $INPUT | tail -1| cut -f2 -d"="`
  GM_evGW=`grep "GM@evGW  correlation energy =" $INPUT | tail -1 | cut -f2 -d"="`
  BSE1_evGW=`grep "BSE@evGW correlation energy (singlet)" $INPUT | cut -f2 -d"="`
  BSE3_evGW=`grep "BSE@evGW correlation energy (triplet)" $INPUT | cut -f2 -d"="`

  echo -e "\t" $RPA_evGW "\t" $GM_evGW "\t" $BSE1_evGW "\t" $BSE3_evGW

  echo 
  echo '*** CIS and TDHF excitation energy (singlet & triplet) ***'
  CIS1=`grep "|     1 |" $INPUT | head -1 | cut -f4 -d"|"`
  CIS3=`grep "|     1 |" $INPUT | head -2 | cut -f4 -d"|" | tail -1`
  TDHF1=`grep "|     1 |" $INPUT | head -3 | cut -f4 -d"|" | tail -1`
  TDHF3=`grep "|     1 |" $INPUT | head -4 | cut -f4 -d"|" | tail -1`
  echo -e "\t" $CIS1 "\t" $CIS3 "\t" $TDHF1 "\t" $TDHF3

  echo 
  echo '*** BSE@G0W0 and BSE@evGW excitation energy (singlet & triplet) ***'
  G0W01=`grep "|     1 |" $INPUT | head -6 | cut -f4 -d"|" | tail -1`
  G0W03=`grep "|     1 |" $INPUT | head -7 | cut -f4 -d"|" | tail -1`
  evGW1=`grep "|     1 |" $INPUT | tail -2 | head -1 | cut -f4 -d"|"`
  evGW3=`grep "|     1 |" $INPUT | tail -1 | cut -f4 -d"|"`
  echo -e "\t" $G0W01 "\t" $G0W03 "\t" $evGW1 "\t" $evGW3

  echo 
  echo '*** DONE ***'

