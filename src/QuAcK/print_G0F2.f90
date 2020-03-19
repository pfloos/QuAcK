subroutine print_G0F2(nBas,nO,eHF,eGF2)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: eGF2(nBas)

  integer                            :: x
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGF2(LUMO)-eGF2(HOMO)

! Dump results

  write(*,*)'-------------------------------------------'
  write(*,*)' Frequency-dependent G0F2 calculation'
  write(*,*)'-------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','e_G0F2 (eV)','|'
  write(*,*)'-------------------------------------------'

  do x=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',x,'|',eHF(x)*HaToeV,'|',eGF2(x)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------'
  write(*,'(2X,A27,F15.6)') 'G0F2  HOMO      energy (eV):',eGF2(HOMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'G0F2  LUMO      energy (eV):',eGF2(LUMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'G0F2  HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------'
  write(*,*)

end subroutine print_G0F2
