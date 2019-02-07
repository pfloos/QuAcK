subroutine ReadGeometry(NAtoms,ZNuc,XYZAtoms)

! Read molecular geometry

  implicit none

! Input variables
  integer,intent(in)            :: NAtoms
  double precision,intent(out)  :: ZNuc(NAtoms),XYZAtoms(NAtoms,3)

! Local variables
  integer                       :: i

! Open file with geometry specification
  open(unit=1,file='input/molecule')

! Read number of atoms
  read(1,*)
  read(1,*)
  read(1,*)

  do i=1,NAtoms
    read(1,*) ZNuc(i),XYZAtoms(i,1),XYZAtoms(i,2),XYZAtoms(i,3)
  enddo

! Print geometry
  write(*,'(A28)') 'Molecular geometry'
  write(*,'(A28)') '------------------'
  do i=1,NAtoms
    write(*,'(A28,1X,I16)') 'Atom n. ',i
    write(*,'(A28,1X,F16.10)') 'Z = ',ZNuc(i)
    write(*,'(A28,1X,F16.10,F16.10,F16.10)') 'Atom coordinates:',XYZAtoms(i,1),XYZAtoms(i,2),XYZAtoms(i,3)
  enddo
  write(*,'(A28)') '------------------'
  write(*,*)

! Close file with geometry specification
  close(unit=1)

end subroutine ReadGeometry
