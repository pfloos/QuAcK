subroutine read_geometry(working_dir,nNuc,ZNuc,rNuc,ENuc)

! Read molecular geometry

  implicit none
  include 'parameters.h'

! Input variables

  character(len=256),intent(in) :: working_dir

! Ouput variables

  integer,intent(in)            :: nNuc

! Local variables

  integer                       :: i,j
  double precision              :: RAB
  character(len=2)              :: El
  integer,external              :: element_number

! Ouput variables

  double precision,intent(out)  :: ZNuc(nNuc),rNuc(nNuc,ncart),ENuc

  integer                       :: status
  character(len=256)            :: file_path

! Open file with geometry specification

  file_path = trim(working_dir) // '/input/molecule'
  open(unit=10, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

      ! Read geometry and create xyz file for integrals
      open(unit=11,file=trim(working_dir) // '/input/molecule.xyz')
      
      read(10,*) 
      read(10,*) 
      read(10,*) 
      
      write(11,'(I3)') nNuc
      write(11,*) 
      
      do i=1,nNuc
        read(10,*) El,rNuc(i,1),rNuc(i,2),rNuc(i,3)
        write(11,'(A3,1X,3F16.10)') El,rNuc(i,1)*BoToAn,rNuc(i,2)*BoToAn,rNuc(i,3)*BoToAn
        ZNuc(i) = dble(element_number(El))
      end do

      close(unit=11)

    endif

  close(unit=10)

  ! ---

  file_path = trim(working_dir) // '/int/ENuc.dat'
  open(unit=3, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

      read(3,*) ENuc

    endif

  close(unit=3)
      
  ! Compute nuclear repulsion energy
  !ENuc = 0.d0
  !do i=1,nNuc-1
  !  do j=i+1,nNuc
  !    RAB = (rNuc(i,1)-rNuc(j,1))**2 + (rNuc(i,2)-rNuc(j,2))**2 + (rNuc(i,3)-rNuc(j,3))**2
  !    ENuc = ENuc + ZNuc(i)*ZNuc(j)/(AntoBo*sqrt(RAB))
  !  end do
  !end do


! Print geometry
  write(*,'(A28)') '------------------'
  write(*,'(A28)') 'Molecular geometry'
  write(*,'(A28)') '------------------'
  do i=1,nNuc
    write(*,'(A28,1X,I16)') 'Atom n. ',i
    write(*,'(A28,1X,F16.10)') 'Z = ',ZNuc(i)
    write(*,'(A28,1X,F16.10,F16.10,F16.10)') 'Atom coordinates:',(rNuc(i,j),j=1,ncart)
  end do
  write(*,*)
  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,F16.10)') 'Nuclear repulsion energy = ',ENuc
  write(*,'(A28)') '------------------'
  write(*,*)

end subroutine 
