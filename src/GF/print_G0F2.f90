subroutine print_G0F2(nBas,nO,eHF,Sig,eGF2,Z,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: Sig(nBas)
  double precision,intent(in)        :: eGF2(nBas)
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
  Gap = eGF2(LUMO) - eGF2(HOMO)

! Dump results

  write(*,*)'--------------------------------------------------------------------------'
  write(*,*)' Frequency-dependent G0F2 calculation'
  write(*,*)'--------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma (eV)','|','e_G0F2 (eV)','|','Z','|'
  write(*,*)'--------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Sig(p)*HaToeV,'|',eGF2(p)*HaToeV,'|',Z(p),'|'
  enddo

  write(*,*)'--------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'G0F2  HOMO      energy (eV):',eGF2(HOMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0F2  LUMO      energy (eV):',eGF2(LUMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0F2  HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'--------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'G0F2 total energy      :',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A30,F15.6,A3)') 'G0F2 correlation energy:',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
