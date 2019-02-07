subroutine ReadNAtoms(NAtoms)

! Read number of atoms

  implicit none

! Input variables
  integer,intent(out)           :: NAtoms

! Open file with geometry specification
  open(unit=1,file='input/molecule')

! Read number of atoms
  read(1,*)
  read(1,*) NAtoms

! Close file with geometry specification
  close(unit=1)

end subroutine ReadNAtoms
