subroutine print_evRGF2(nBas,nOrb,nC,nO,nV,nR,nSCF,Conv,eHF,SigC,Z,eGF,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nSCF
  double precision,intent(in)   :: Conv
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: SigC(nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Z(nOrb)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: Ec

  integer                       :: i,a
  integer                       :: HOMO
  integer                       :: LUMO
  double precision              :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGF(LUMO) - eGF(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Self-consistent evGF2 calculation'
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
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evRGF2 HOMO      energy =',eGF(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evRGF2 LUMO      energy =',eGF(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evRGF2 HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evRGF2 total energy       =',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A60,F15.6,A3)') 'evRGF2 correlation energy =',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)


end subroutine 
