#! /bin/bash

INPUT=$1

  echo 
  echo '******************************************'
  echo '*** Extracting information of' $INPUT ' ***'
  echo '******************************************'
  echo 

  echo 
  echo '*** WFT information ***'
  grep "Hartree-Fock energy" $INPUT
  EHF=`grep "Hartree-Fock energy" $INPUT | cut -f2 -d":"`
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
  echo '*** RPA@TDHF information: RPA, RPA1 & RPA3  ***'
  RPA_TDHF=`grep "RPA@TDHF correlation energy           =" $INPUT| cut -f2 -d"="`
  RPA1_TDHF=`grep "RPA@TDHF correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  RPA3_TDHF=`grep "RPA@TDHF correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $RPA_TDHF "\t" $RPA1_TDHF "\t" $RPA3_TDHF

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
  CIS1_1=`grep "|     1 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS3_1=`grep "|     1 |" $INPUT | head -2 | cut -f3 -d"|" | tail -1`
  TDHF1_1=`grep "|     1 |" $INPUT | head -3 | cut -f3 -d"|" | tail -1`
  TDHF3_1=`grep "|     1 |" $INPUT | head -4 | cut -f3 -d"|" | tail -1`
  echo -e "\t" $CIS1_1 "\t" $CIS3_1 "\t" $TDHF1_1 "\t" $TDHF3_1

  CIS1_2=`grep "|     2 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS3_2=`grep "|     2 |" $INPUT | head -2 | cut -f3 -d"|" | tail -1`
  TDHF1_2=`grep "|     2 |" $INPUT | head -3 | cut -f3 -d"|" | tail -1`
  TDHF3_2=`grep "|     2 |" $INPUT | head -4 | cut -f3 -d"|" | tail -1`
  echo -e "\t" $CIS1_2 "\t" $CIS3_2 "\t" $TDHF1_2 "\t" $TDHF3_2

  CIS1_3=`grep "|     3 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS3_3=`grep "|     3 |" $INPUT | head -2 | cut -f3 -d"|" | tail -1`
  TDHF1_3=`grep "|     3 |" $INPUT | head -3 | cut -f3 -d"|" | tail -1`
  TDHF3_3=`grep "|     3 |" $INPUT | head -4 | cut -f3 -d"|" | tail -1`
  echo -e "\t" $CIS1_3 "\t" $CIS3_3 "\t" $TDHF1_3 "\t" $TDHF3_3

  CIS1_4=`grep "|     4 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS3_4=`grep "|     4 |" $INPUT | head -2 | cut -f3 -d"|" | tail -1`
  TDHF1_4=`grep "|     4 |" $INPUT | head -3 | cut -f3 -d"|" | tail -1`
  TDHF3_4=`grep "|     4 |" $INPUT | head -4 | cut -f3 -d"|" | tail -1`
  echo -e "\t" $CIS1_4 "\t" $CIS3_4 "\t" $TDHF1_4 "\t" $TDHF3_4

  CIS1_5=`grep "|     5 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS3_5=`grep "|     5 |" $INPUT | head -2 | cut -f3 -d"|" | tail -1`
  TDHF1_5=`grep "|     5 |" $INPUT | head -3 | cut -f3 -d"|" | tail -1`
  TDHF3_5=`grep "|     5 |" $INPUT | head -4 | cut -f3 -d"|" | tail -1`
  echo -e "\t" $CIS1_5 "\t" $CIS3_5 "\t" $TDHF1_5 "\t" $TDHF3_5

  echo 
  echo '*** BSE@G0W0 and BSE@evGW excitation energy (singlet & triplet) ***'
  G0W01_1=`grep "|     1 |" $INPUT | head -6 | cut -f3 -d"|" | tail -1`
  G0W03_1=`grep "|     1 |" $INPUT | head -7 | cut -f3 -d"|" | tail -1`
  evGW1_1=`grep "|     1 |" $INPUT | tail -2 | head -1 | cut -f3 -d"|"`
  evGW3_1=`grep "|     1 |" $INPUT | tail -1 | cut -f3 -d"|"`
  echo -e "\t" $G0W01_1 "\t" $G0W03_1 "\t" $evGW1_1 "\t" $evGW3_1

  G0W01_2=`grep "|     2 |" $INPUT | head -6 | cut -f3 -d"|" | tail -1`
  G0W03_2=`grep "|     2 |" $INPUT | head -7 | cut -f3 -d"|" | tail -1`
  evGW1_2=`grep "|     2 |" $INPUT | tail -2 | head -1 | cut -f3 -d"|"`
  evGW3_2=`grep "|     2 |" $INPUT | tail -1 | cut -f3 -d"|"`
  echo -e "\t" $G0W01_2 "\t" $G0W03_2 "\t" $evGW1_2 "\t" $evGW3_2

  G0W01_3=`grep "|     3 |" $INPUT | head -6 | cut -f3 -d"|" | tail -1`
  G0W03_3=`grep "|     3 |" $INPUT | head -7 | cut -f3 -d"|" | tail -1`
  evGW1_3=`grep "|     3 |" $INPUT | tail -2 | head -1 | cut -f3 -d"|"`
  evGW3_3=`grep "|     3 |" $INPUT | tail -1 | cut -f3 -d"|"`
  echo -e "\t" $G0W01_3 "\t" $G0W03_3 "\t" $evGW1_3 "\t" $evGW3_3

  G0W01_4=`grep "|     4 |" $INPUT | head -6 | cut -f3 -d"|" | tail -1`
  G0W03_4=`grep "|     4 |" $INPUT | head -7 | cut -f3 -d"|" | tail -1`
  evGW1_4=`grep "|     4 |" $INPUT | tail -2 | head -1 | cut -f3 -d"|"`
  evGW3_4=`grep "|     4 |" $INPUT | tail -1 | cut -f3 -d"|"`
  echo -e "\t" $G0W01_4 "\t" $G0W03_4 "\t" $evGW1_4 "\t" $evGW3_4

  G0W01_5=`grep "|     5 |" $INPUT | head -6 | cut -f3 -d"|" | tail -1`
  G0W03_5=`grep "|     5 |" $INPUT | head -7 | cut -f3 -d"|" | tail -1`
  evGW1_5=`grep "|     5 |" $INPUT | tail -2 | head -1 | cut -f3 -d"|"`
  evGW3_5=`grep "|     5 |" $INPUT | tail -1 | cut -f3 -d"|"`
  echo -e "\t" $G0W01_5 "\t" $G0W03_5 "\t" $evGW1_5 "\t" $evGW3_5

  echo '*** MATHEMATICA OUTPUT ***'
  echo -e "\t" $EHF "\t" $CIS1_1 "\t" $CIS1_2 "\t" $CIS1_3 "\t" $CIS1_4 "\t" $CIS1_5 "\t" $CIS3_1 "\t" $CIS3_2 "\t" $CIS3_3 "\t" $CIS3_4 "\t" $CIS3_5 "\t" $RPA_TDHF "\t" $RPA1_TDHF "\t" $RPA3_TDHF "\t" $TDHF1_1 "\t" $TDHF1_2 "\t" $TDHF1_3 "\t" $TDHF1_4 "\t" $TDHF1_5 "\t" $TDHF3_1 "\t" $TDHF3_2 "\t" $TDHF3_3 "\t" $TDHF3_4 "\t" $TDHF3_5 "\t" $RPA_G0W0 "\t" $GM_G0W0 "\t" $BSE1_G0W0 "\t" $BSE3_G0W0 "\t" $G0W01_1 "\t" $G0W01_2 "\t" $G0W01_3 "\t" $G0W01_4 "\t" $G0W01_5 "\t" $G0W03_1 "\t" $G0W03_2 "\t" $G0W03_3 "\t" $G0W03_4 "\t" $G0W03_5 "\t" $RPA_evGW "\t" $GM_evGW "\t" $BSE1_evGW "\t" $BSE3_evGW "\t" $evGW1_1 "\t" $evGW1_2 "\t" $evGW1_3 "\t" $evGW1_4 "\t" $evGW1_5 "\t" $evGW3_1 "\t" $evGW3_2 "\t" $evGW3_3 "\t" $evGW3_4 "\t" $evGW3_5 
  echo '*** DONE ***'

