subroutine print_RG0F2(nOrb,nC,nO,nV,nR,eHF,SigC,eGF,Z,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: SigC(nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Z(nOrb)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: Ec

  integer                       :: i,a
  double precision              :: eHOMO,eLUMO,Gap

! HOMO and LUMO

  eHOMO = maxval(eGF(1:nO))
  eLUMO = minval(eGF(nO+1:nOrb))
  Gap   = eLUMO - eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' G0F2@RHF calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GF2 (eV)','|','Z','|','e_GF2 (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',i,'|',eHF(i)*HaToeV,'|',SigC(i)*HaToeV,'|',Z(i),'|',eGF(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',a,'|',eHF(a)*HaToeV,'|',SigC(a)*HaToeV,'|',Z(a),'|',eGF(a)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  HOMO      energy =',eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  LUMO      energy =',eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2 total energy       =',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2 correlation energy =',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
