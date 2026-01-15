subroutine print_MOMRHF(nBas, nOrb, nO, eHF, cHF, ENuc, ET, EV, EJ, EK, ERHF, dipole, occupations)

! Print one-electron energies and other stuff for MOMRHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas, nOrb
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nOrb)
  double precision,intent(in)        :: cHF(nBas,nOrb)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: dipole(ncart)
  integer,intent(in)                 :: occupations(nO)

! Local variables

  integer                            :: ixyz
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: S,S2
  integer,allocatable                :: virtual(:)

  logical                            :: dump_orb = .false.

  allocate(virtual(nOrb-nO))

! HOMO and LUMO

  call non_occupied(nO, nOrb, occupations,virtual)  
  HOMO = maxval(occupations)
  LUMO = minval(virtual)
  Gap = eHF(LUMO)-eHF(HOMO)

  S2 = 0d0
  S  = 0d0

! Dump results

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',ERHF,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' RHF          energy = ',ERHF + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6,A3)')  ' HF HOMO      energy = ',eHF(HOMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' HF LUMO      energy = ',eHF(LUMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' HF HOMO-LUMO gap    = ',Gap*HaToeV,' eV'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6)')     ' <Sz>                = ',S
  write(*,'(A33,1X,F16.6)')     ' <S^2>               = ',S2
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36)')           ' Dipole moment (Debye)    '
  write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
  write(*,'(10X,4F10.4)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A50)')           '---------------------------------------'
  write(*,*)

! Print results

  if(dump_orb) then 
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A50)') ' RHF orbital coefficients '
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas, nOrb, cHF)
    write(*,*)
  end if

  write(*,'(A50)') '-------------------------------------------'
  write(*,'(A50)') ' Doubly occupied RHF orbital energies (au) '
  write(*,'(A50)') '-------------------------------------------'
  call vecout(nO, eHF(occupations(:)))
  write(*,*)
  
  write(*,'(A50)') '-----------------------------------'
  write(*,'(A50)') ' Virtual RHF orbital energies (au) '
  write(*,'(A50)') '-----------------------------------'
  call vecout(nOrb - nO, eHF(virtual(:)))
  write(*,*)
  
  print *, "Orbital occupations for MOMRHF:"
  print *, occupations(:)
  print *, ""
  
  deallocate(virtual)

end subroutine 
