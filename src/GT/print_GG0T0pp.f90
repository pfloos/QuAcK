subroutine print_GG0T0pp(nOrb,nC,nO,nV,nR,eHF,ENuc,EGHF,SigC,Z,eGT,EcGM,EcRPA)

! Print one-electron energies and other stuff for G0T0

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: EcGM
  double precision,intent(in)   :: EcRPA
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: SigC(nOrb)
  double precision,intent(in)   :: Z(nOrb)
  double precision,intent(in)   :: eGT(nOrb)

  integer                       :: i,a
  double precision              :: eHOMO,eLUMO,Gap

! HOMO and LUMO

  eHOMO = maxval(eGT(1:nO))
  eLUMO = minval(eGT(nO+1:nOrb))
  Gap   = eLUMO - eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' G0T0pp@GHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GTpp (eV)','|','Z','|','e_GTpp (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',i,'|',eHF(i)*HaToeV,'|',SigC(i)*HaToeV,'|',Z(i),'|',eGT(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',a,'|',eHF(a)*HaToeV,'|',SigC(a)*HaToeV,'|',Z(a),'|',eGT(a)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@GHF HOMO      energy                 = ',eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@GHF LUMO      energy                 = ',eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@GHF HOMO-LUMO gap                    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@G0T0pp@GHF correlation energy           = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@G0T0pp@GHF total       energy           = ',ENuc + EGHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0T0pp@GHF correlation energy           = ',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0T0pp@GHF total       energy           = ',ENuc + EGHF + EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
