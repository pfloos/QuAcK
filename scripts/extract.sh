#! /bin/bash

INPUT=$1

  echo 
  echo '******************************************'
  echo '*** Extracting information of' $INPUT ' ***'
  echo '******************************************'
  echo 

  echo 
  echo '*** WFT information ***'
  echo 
  grep "Hartree-Fock energy" $INPUT
  EHF=`grep "Hartree-Fock energy" $INPUT | cut -f2 -d"="`
  grep "MP2 correlation energy" $INPUT
  EcMP2=`grep "MP2 correlation energy" $INPUT | cut -f2 -d"="`
  grep "Ec(MP2) =" $INPUT
  grep "Ec(CCD) =" $INPUT
  grep "Ec(CCSD) =" $INPUT
  grep "Ec(CCSD(T)) =" $INPUT

#  echo 
#  echo '*** Gap information: HF, G0F2, GF2, G0W0 & evGW ***'
#  HF=`grep "HF HOMO-LUMO gap    (eV):" $INPUT | cut -f2 -d":"`
#  G0F2=`grep "GF2  HOMO-LUMO gap    (eV):" $INPUT | head -1 | cut -f2 -d":"`
#  GF2=`grep "GF2  HOMO-LUMO gap    (eV):" $INPUT | tail -1 | cut -f2 -d":"`
#  G0W0=`grep "G0W0 HOMO-LUMO gap    (eV):" $INPUT | cut -f2 -d":"`
#  evGW=`grep "evGW HOMO-LUMO gap    (eV):" $INPUT | tail -1 | cut -f2 -d":"`

#  echo -e "\t" $HF "\t" $G0F2 "\t" $GF2 "\t" $G0W0 "\t" $evGW

  echo 
  echo '*** RPA information: Tr@RPA (singlet), Tr@RPA (triplet), AC@RPA (singlet), AC@RPA (triplet) ***'
  echo 
  Tr_RPA_1=`grep "Tr@RPA  correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_RPA_3=`grep "Tr@RPA  correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  AC_RPA_1=`grep "AC@RPA  correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  AC_RPA_3=`grep "AC@RPA  correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $Tr_RPA_1 "\t" $Tr_RPA_3 "\t" $AC_RPA_1 "\t" $AC_RPA_3

  echo 
  echo '*** RPAx information: Tr@RPAx (singlet), Tr@RPAx (triplet), AC@RPAx (singlet), AC@RPAx (triplet) ***'
  echo 
  Tr_RPAx_1=`grep "Tr@RPAx correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_RPAx_3=`grep "Tr@RPAx correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  AC_RPAx_1=`grep "AC@RPAx correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  AC_RPAx_3=`grep "AC@RPAx correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $Tr_RPAx_1 "\t" $Tr_RPAx_3 "\t" $AC_RPAx_1 "\t" $AC_RPAx_3

  echo 
  echo '*** G0W0 information: Tr@RPA (singlet), Tr@RPA (triplet), Tr@BSE (singlet), Tr@BSE (triplet), AC@BSE (singlet), AC@BSE (triplet) ***'
  echo 
  Tr_RPA_G0W0_1=`grep "Tr@RPA@G0W0 correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_RPA_G0W0_3=`grep "Tr@RPA@G0W0 correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_G0W0_1=`grep "Tr@BSE@G0W0 correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_G0W0_3=`grep "Tr@BSE@G0W0 correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_G0W0_1=`grep "AC@BSE@G0W0 correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_G0W0_3=`grep "AC@BSE@G0W0 correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $Tr_RPA_G0W0_1 "\t" $Tr_RPA_G0W0_3 "\t" $Tr_BSE_G0W0_1 "\t" $Tr_BSE_G0W0_3 "\t" $AC_BSE_G0W0_1 "\t" $AC_BSE_G0W0_3

  echo 
  echo '*** evGW information: Tr@RPA (singlet), Tr@RPA (triplet), Tr@BSE (singlet), Tr@BSE (triplet), AC@BSE (singlet), AC@BSE (triplet) ***'
  echo 
  Tr_RPA_evGW_1=`grep "Tr@RPA@evGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_RPA_evGW_3=`grep "Tr@RPA@evGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_evGW_1=`grep "Tr@BSE@evGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_evGW_3=`grep "Tr@BSE@evGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_evGW_1=`grep "AC@BSE@evGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_evGW_3=`grep "AC@BSE@evGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $Tr_RPA_evGW_1 "\t" $Tr_RPA_evGW_3 "\t" $Tr_BSE_evGW_1 "\t" $Tr_BSE_evGW_3 "\t" $AC_BSE_evGW_1 "\t" $AC_BSE_evGW_3


  echo 
  echo '*** qsGW information: Tr@RPA (singlet), Tr@RPA (triplet), Tr@BSE (singlet), Tr@BSE (triplet), AC@BSE (singlet), AC@BSE (triplet) ***'
  echo 
  Tr_RPA_qsGW_1=`grep "Tr@RPA@qsGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_RPA_qsGW_3=`grep "Tr@RPA@qsGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_qsGW_1=`grep "Tr@BSE@qsGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  Tr_BSE_qsGW_3=`grep "Tr@BSE@qsGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_qsGW_1=`grep "AC@BSE@qsGW correlation energy (singlet) =" $INPUT| cut -f2 -d"="`
  AC_BSE_qsGW_3=`grep "AC@BSE@qsGW correlation energy (triplet) =" $INPUT| cut -f2 -d"="`

  echo -e "\t" $Tr_RPA_qsGW_1 "\t" $Tr_RPA_qsGW_3 "\t" $Tr_BSE_qsGW_1 "\t" $Tr_BSE_qsGW_3 "\t" $AC_BSE_qsGW_1 "\t" $AC_BSE_qsGW_3

  echo 
  echo '*** CIS excitation energy (singlet & triplet) ***'
  echo 

  CIS_1_1=`grep "|     1 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS_1_2=`grep "|     2 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS_1_3=`grep "|     3 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS_1_4=`grep "|     4 |" $INPUT | head -1 | cut -f3 -d"|"`
  CIS_1_5=`grep "|     5 |" $INPUT | head -1 | cut -f3 -d"|"`

  CIS_3_1=`grep "|     1 |" $INPUT | head -2 | tail -1 | cut -f3 -d"|"`
  CIS_3_2=`grep "|     2 |" $INPUT | head -2 | tail -1 | cut -f3 -d"|"`
  CIS_3_3=`grep "|     3 |" $INPUT | head -2 | tail -1 | cut -f3 -d"|"`
  CIS_3_4=`grep "|     4 |" $INPUT | head -2 | tail -1 | cut -f3 -d"|"`
  CIS_3_5=`grep "|     5 |" $INPUT | head -2 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $CIS_1_1 "\t" $CIS_3_1 
  echo -e "\t" $CIS_1_2 "\t" $CIS_3_2
  echo -e "\t" $CIS_1_3 "\t" $CIS_3_3
  echo -e "\t" $CIS_1_4 "\t" $CIS_3_4
  echo -e "\t" $CIS_1_5 "\t" $CIS_3_5

  echo 
  echo '*** RPA excitation energy (singlet & triplet) ***'
  echo 

  RPA_1_1=`grep "|     1 |" $INPUT | head -3 | tail -1 | cut -f3 -d"|"`
  RPA_1_2=`grep "|     2 |" $INPUT | head -3 | tail -1 | cut -f3 -d"|"`
  RPA_1_3=`grep "|     3 |" $INPUT | head -3 | tail -1 | cut -f3 -d"|"`
  RPA_1_4=`grep "|     4 |" $INPUT | head -3 | tail -1 | cut -f3 -d"|"`
  RPA_1_5=`grep "|     5 |" $INPUT | head -3 | tail -1 | cut -f3 -d"|"`

  RPA_3_1=`grep "|     1 |" $INPUT | head -4 | tail -1 | cut -f3 -d"|"`
  RPA_3_2=`grep "|     2 |" $INPUT | head -4 | tail -1 | cut -f3 -d"|"`
  RPA_3_3=`grep "|     3 |" $INPUT | head -4 | tail -1 | cut -f3 -d"|"`
  RPA_3_4=`grep "|     4 |" $INPUT | head -4 | tail -1 | cut -f3 -d"|"`
  RPA_3_5=`grep "|     5 |" $INPUT | head -4 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $RPA_1_1 "\t" $RPA_3_1 
  echo -e "\t" $RPA_1_2 "\t" $RPA_3_2
  echo -e "\t" $RPA_1_3 "\t" $RPA_3_3
  echo -e "\t" $RPA_1_4 "\t" $RPA_3_4
  echo -e "\t" $RPA_1_5 "\t" $RPA_3_5

  echo 
  echo '*** RPAx excitation energy (singlet & triplet) ***'
  echo 

  RPAx_1_1=`grep "|     1 |" $INPUT | head -5 | tail -1 | cut -f3 -d"|"`
  RPAx_1_2=`grep "|     2 |" $INPUT | head -5 | tail -1 | cut -f3 -d"|"`
  RPAx_1_3=`grep "|     3 |" $INPUT | head -5 | tail -1 | cut -f3 -d"|"`
  RPAx_1_4=`grep "|     4 |" $INPUT | head -5 | tail -1 | cut -f3 -d"|"`
  RPAx_1_5=`grep "|     5 |" $INPUT | head -5 | tail -1 | cut -f3 -d"|"`

  RPAx_3_1=`grep "|     1 |" $INPUT | head -6 | tail -1 | cut -f3 -d"|"`
  RPAx_3_2=`grep "|     2 |" $INPUT | head -6 | tail -1 | cut -f3 -d"|"`
  RPAx_3_3=`grep "|     3 |" $INPUT | head -6 | tail -1 | cut -f3 -d"|"`
  RPAx_3_4=`grep "|     4 |" $INPUT | head -6 | tail -1 | cut -f3 -d"|"`
  RPAx_3_5=`grep "|     5 |" $INPUT | head -6 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $RPAx_1_1 "\t" $RPAx_3_1 
  echo -e "\t" $RPAx_1_2 "\t" $RPAx_3_2
  echo -e "\t" $RPAx_1_3 "\t" $RPAx_3_3
  echo -e "\t" $RPAx_1_4 "\t" $RPAx_3_4
  echo -e "\t" $RPAx_1_5 "\t" $RPAx_3_5

  echo 
  echo '*** BSE@G0W0 excitation energy (singlet & triplet) ***'
  echo 

  G0W0_1_1=`grep "|     1 |" $INPUT | head -7 | tail -1 | cut -f3 -d"|"`
  G0W0_1_2=`grep "|     2 |" $INPUT | head -7 | tail -1 | cut -f3 -d"|"`
  G0W0_1_3=`grep "|     3 |" $INPUT | head -7 | tail -1 | cut -f3 -d"|"`
  G0W0_1_4=`grep "|     4 |" $INPUT | head -7 | tail -1 | cut -f3 -d"|"`
  G0W0_1_5=`grep "|     5 |" $INPUT | head -7 | tail -1 | cut -f3 -d"|"`

  G0W0_3_1=`grep "|     1 |" $INPUT | head -8 | tail -1 | cut -f3 -d"|"`
  G0W0_3_2=`grep "|     2 |" $INPUT | head -8 | tail -1 | cut -f3 -d"|"`
  G0W0_3_3=`grep "|     3 |" $INPUT | head -8 | tail -1 | cut -f3 -d"|"`
  G0W0_3_4=`grep "|     4 |" $INPUT | head -8 | tail -1 | cut -f3 -d"|"`
  G0W0_3_5=`grep "|     5 |" $INPUT | head -8 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $G0W0_1_1 "\t" $G0W0_3_1 
  echo -e "\t" $G0W0_1_2 "\t" $G0W0_3_2
  echo -e "\t" $G0W0_1_3 "\t" $G0W0_3_3
  echo -e "\t" $G0W0_1_4 "\t" $G0W0_3_4
  echo -e "\t" $G0W0_1_5 "\t" $G0W0_3_5

  echo 
  echo '*** BSE@evGW excitation energy (singlet & triplet) ***'
  echo 

  evGW_1_1=`grep "|     1 |" $INPUT | head -9 | tail -1 | cut -f3 -d"|"`
  evGW_1_2=`grep "|     2 |" $INPUT | head -9 | tail -1 | cut -f3 -d"|"`
  evGW_1_3=`grep "|     3 |" $INPUT | head -9 | tail -1 | cut -f3 -d"|"`
  evGW_1_4=`grep "|     4 |" $INPUT | head -9 | tail -1 | cut -f3 -d"|"`
  evGW_1_5=`grep "|     5 |" $INPUT | head -9 | tail -1 | cut -f3 -d"|"`

  evGW_3_1=`grep "|     1 |" $INPUT | head -10 | tail -1 | cut -f3 -d"|"`
  evGW_3_2=`grep "|     2 |" $INPUT | head -10 | tail -1 | cut -f3 -d"|"`
  evGW_3_3=`grep "|     3 |" $INPUT | head -10 | tail -1 | cut -f3 -d"|"`
  evGW_3_4=`grep "|     4 |" $INPUT | head -10 | tail -1 | cut -f3 -d"|"`
  evGW_3_5=`grep "|     5 |" $INPUT | head -10 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $evGW_1_1 "\t" $evGW_3_1 
  echo -e "\t" $evGW_1_2 "\t" $evGW_3_2
  echo -e "\t" $evGW_1_3 "\t" $evGW_3_3
  echo -e "\t" $evGW_1_4 "\t" $evGW_3_4
  echo -e "\t" $evGW_1_5 "\t" $evGW_3_5

  echo 
  echo '*** BSE@qsGW excitation energy (singlet & triplet) ***'
  echo 

  qsGW_1_1=`grep "|     1 |" $INPUT | head -11 | tail -1 | cut -f3 -d"|"`
  qsGW_1_2=`grep "|     2 |" $INPUT | head -11 | tail -1 | cut -f3 -d"|"`
  qsGW_1_3=`grep "|     3 |" $INPUT | head -11 | tail -1 | cut -f3 -d"|"`
  qsGW_1_4=`grep "|     4 |" $INPUT | head -11 | tail -1 | cut -f3 -d"|"`
  qsGW_1_5=`grep "|     5 |" $INPUT | head -11 | tail -1 | cut -f3 -d"|"`

  qsGW_3_1=`grep "|     1 |" $INPUT | head -12 | tail -1 | cut -f3 -d"|"`
  qsGW_3_2=`grep "|     2 |" $INPUT | head -12 | tail -1 | cut -f3 -d"|"`
  qsGW_3_3=`grep "|     3 |" $INPUT | head -12 | tail -1 | cut -f3 -d"|"`
  qsGW_3_4=`grep "|     4 |" $INPUT | head -12 | tail -1 | cut -f3 -d"|"`
  qsGW_3_5=`grep "|     5 |" $INPUT | head -12 | tail -1 | cut -f3 -d"|"`

  echo -e "\t" $qsGW_1_1 "\t" $qsGW_3_1 
  echo -e "\t" $qsGW_1_2 "\t" $qsGW_3_2
  echo -e "\t" $qsGW_1_3 "\t" $qsGW_3_3
  echo -e "\t" $qsGW_1_4 "\t" $qsGW_3_4
  echo -e "\t" $qsGW_1_5 "\t" $qsGW_3_5

  echo 
  echo '*** MATHEMATICA OUTPUT ***'
  echo 
  echo -e "\t" $EHF "\t" $EcMP2 "\t" $Tr_RPA_1 "\t" $Tr_RPA_3 "\t" $AC_RPA_1 "\t" $AC_RPA_3 "\t" $Tr_RPAx_1 "\t" $Tr_RPAx_3 "\t" $AC_RPAx_1 "\t" $AC_RPAx_3 "\t" $Tr_RPA_G0W0_1 "\t" $Tr_RPA_G0W0_3 "\t" $Tr_BSE_G0W0_1 "\t" $Tr_BSE_G0W0_3 "\t" $AC_BSE_G0W0_1 "\t" $AC_BSE_G0W0_3 "\t" $CIS_1_1 "\t" $CIS_1_2 "\t" $CIS_1_3 "\t" $CIS_1_4 "\t" $CIS_1_5 "\t" $CIS_3_1 "\t" $CIS_3_2 "\t" $CIS_3_3 "\t" $CIS_3_4 "\t" $CIS_3_5 "\t" $RPA_1_1 "\t" $RPA_1_2 "\t" $RPA_1_3 "\t" $RPA_1_4 "\t" $RPA_1_5 "\t" $RPA_3_1 "\t" $RPA_3_2 "\t" $RPA_3_3 "\t" $RPA_3_4 "\t" $RPA_3_5 "\t" $RPAx_1_1 "\t" $RPAx_1_2 "\t" $RPAx_1_3 "\t" $RPAx_1_4 "\t" $RPAx_1_5 "\t" $RPAx_3_1 "\t" $RPAx_3_2 "\t" $RPAx_3_3 "\t" $RPAx_3_4 "\t" $RPAx_3_5 "\t" $G0W0_1_1 "\t" $G0W0_1_2 "\t" $G0W0_1_3 "\t" $G0W0_1_4 "\t" $G0W0_1_5 "\t" $G0W0_3_1 "\t" $G0W0_3_2 "\t" $G0W0_3_3 "\t" $G0W0_3_4 "\t" $G0W0_3_5 "\t" $Tr_RPA_evGW_1 "\t" $Tr_RPA_evGW_3 "\t" $Tr_BSE_evGW_1 "\t" $Tr_BSE_evGW_3 "\t" $AC_BSE_evGW_1 "\t" $AC_BSE_evGW_3 "\t" $evGW_1_1 "\t" $evGW_1_2 "\t" $evGW_1_3 "\t" $evGW_1_4 "\t" $evGW_1_5 "\t" $evGW_3_1 "\t" $evGW_3_2 "\t" $evGW_3_3 "\t" $evGW_3_4 "\t" $evGW_3_5 "\t" $Tr_RPA_qsGW_1 "\t" $Tr_RPA_qsGW_3 "\t" $Tr_BSE_qsGW_1 "\t" $Tr_BSE_qsGW_3 "\t" $AC_BSE_qsGW_1 "\t" $AC_BSE_qsGW_3 "\t" $qsGW_1_1 "\t" $qsGW_1_2 "\t" $qsGW_1_3 "\t" $qsGW_1_4 "\t" $qsGW_1_5 "\t" $qsGW_3_1 "\t" $qsGW_3_2 "\t" $qsGW_3_3 "\t" $qsGW_3_4 "\t" $qsGW_3_5
  echo 
  echo '*** DONE ***'
  echo 

