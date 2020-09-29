subroutine print_G0F3(nBas,nO,eHF,Z,eGF3)

! Print one-electron energies and other stuff for GF3

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: eHF(nBas),eGF3(nBas),Z(nBas)

  integer                            :: x,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGF3(LUMO)-eGF3(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------'
  write(*,*)' Frequency-dependent G0F3 calculation'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Z','|','e_G0F3 (eV)','|'
  write(*,*)'-------------------------------------------------------------'

  do x=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',x,'|',eHF(x)*HaToeV,'|',Z(x),'|',eGF3(x)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------'
  write(*,'(2X,A27,F15.6)') 'G0F3  HOMO      energy (eV):',eGF3(HOMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'G0F3  LUMO      energy (eV):',eGF3(LUMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'G0F3  HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

end subroutine print_G0F3
