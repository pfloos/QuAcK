
! ---

subroutine print_cRHF(nBas, nOrb, nO, eHF, cHF, ENuc, ET, EV, EW,EJ, EK, ERHF)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas, nOrb
  integer,intent(in)                 :: nO
  complex*16,intent(in)              :: eHF(nOrb)
  complex*16,intent(in)              :: cHF(nBas,nOrb)
  double precision,intent(in)        :: ENuc
  complex*16,intent(in)              :: ET
  complex*16,intent(in)              :: EJ
  complex*16,intent(in)              :: EK
  complex*16,intent(in)              :: EV
  complex*16,intent(in)              :: EW
  complex*16,intent(in)              :: ERHF

! Local variables

  integer                            :: HOMO
  integer                            :: LUMO
  complex*16                         :: Gap
  double precision                   :: S,S2

  logical                            :: dump_orb = .false.

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eHF(LUMO)-eHF(HOMO)

  S2 = 0d0
  S  = 0d0

! Dump results

  write(*,*)
  write(*,'(A69)')           '---------------------------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A69)')           '---------------------------------------------------------'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' One-electron energy = ',real(ET + EV),'+',aimag(ET+EV),'i',' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Kinetic      energy = ',real(ET),'+',aimag(ET),'i',' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Potential    energy = ',real(EV),'+',aimag(EV),'i',' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' CAP          energy = ',real(EW),'+',aimag(EW),'i',' au'
  write(*,'(A69)')           '---------------------------------------------------'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Two-electron energy = ',real(EJ + EK),'+',aimag(EJ+EK),'i',' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Hartree      energy = ',real(EJ),'+',aimag(EJ),'i',' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Exchange     energy = ',real(EK),'+',aimag(EK),'i',' au'
  write(*,'(A69)')           '---------------------------------------------------------'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Electronic   energy = ',real(ERHF),'+',aimag(ERHF),'i',' au'
  write(*,'(A33,1X,F16.10,19X,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,1X,A1,F16.10,A1,A3)') ' cRHF         energy = ',real(ERHF + ENuc),'+',aimag(ERHF+ENuc),'i',' au'
  write(*,'(A69)')           '---------------------------------------------------------'
  write(*,'(A33,1X,F16.6,A3)')  ' HF HOMO      energy = ',real(eHF(HOMO))*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' HF LUMO      energy = ',real(eHF(LUMO))*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' HF HOMO-LUMO gap    = ',real(Gap)*HaToeV,' eV'
  write(*,'(A69)')           '---------------------------------------------------------'
  write(*,'(A33,1X,F16.6)')     ' <Sz>                = ',S
  write(*,'(A33,1X,F16.6)')     ' <S^2>               = ',S2
  write(*,*)

! Print results

  if(dump_orb) then 
    write(*,'(A69)') '---------------------------------------'
    write(*,'(A69)') ' cRHF orbital coefficients '
    write(*,'(A69)') '---------------------------------------'
    call complex_matout(nBas, nOrb, cHF)
    write(*,*)
  end if
  write(*,'(A37)') '-------------------------------'
  write(*,'(A37)') ' cRHF orbital energies (au) '
  write(*,'(A37)') '-------------------------------'
  call complex_vecout(nOrb, eHF)
  write(*,*)

end subroutine 
