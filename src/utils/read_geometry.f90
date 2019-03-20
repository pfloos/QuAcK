subroutine read_geometry(nNuc,ZNuc,rNuc,ENuc)

! Read molecular geometry

  implicit none

  include 'parameters.h'

! Ouput variables

  integer,intent(in)            :: nNuc

! Local variables

  integer                       :: i,j
  double precision              :: RAB
  character(len=2)              :: El
  integer,external              :: element_number

! Ouput variables

  double precision,intent(out)  :: ZNuc(nNuc),rNuc(nNuc,ncart),ENuc

! Open file with geometry specification

  open(unit=1,file='input/molecule')

! Read geometry

  read(1,*) 
  read(1,*) 
  read(1,*) 

  do i=1,nNuc
    read(1,*) El,rNuc(i,1),rNuc(i,2),rNuc(i,3)
    ZNuc(i) = element_number(El)
  enddo

! Compute nuclear repulsion energy

  ENuc = 0

  do i=1,nNuc-1
    do j=i+1,nNuc
      RAB = (rNuc(i,1)-rNuc(j,1))**2 + (rNuc(i,2)-rNuc(j,2))**2 + (rNuc(i,3)-rNuc(j,3))**2
      ENuc = ENuc + ZNuc(i)*ZNuc(j)/sqrt(RAB)
    enddo
  enddo

! Close file with geometry specification
  close(unit=1)

! Print geometry
  write(*,'(A28)') '------------------'
  write(*,'(A28)') 'Molecular geometry'
  write(*,'(A28)') '------------------'
  do i=1,nNuc
    write(*,'(A28,1X,I16)') 'Atom n. ',i
    write(*,'(A28,1X,F16.10)') 'Z = ',ZNuc(i)
    write(*,'(A28,1X,F16.10,F16.10,F16.10)') 'Atom coordinates:',(rNuc(i,j),j=1,ncart)
  enddo
  write(*,*)
  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,F16.10)') 'Nuclear repulsion energy = ',ENuc
  write(*,'(A28)') '------------------'
  write(*,*)

end subroutine read_geometry
