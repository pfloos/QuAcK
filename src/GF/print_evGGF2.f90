subroutine print_evGGF2(nBas,nO,nSCF,Conv,eHF,Sig,Z,eGF,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: Sig(nBas)
  double precision,intent(in)        :: eGF(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: Ec

  integer                            :: p
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap

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

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Sig(p)*HaToeV,'|',Z(p),'|',eGF(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evGGF2 HOMO      energy =',eGF(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGGF2 LUMO      energy =',eGF(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGGF2 HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evGGF2 total energy       =',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A60,F15.6,A3)') 'evGGF2 correlation energy =',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)


end subroutine 
