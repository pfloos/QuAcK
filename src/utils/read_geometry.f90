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

  ! Print geometry

  write(*,*) '============================================================'
  write(*,*) '                    Molecular Geometry                      '
  write(*,*) '============================================================'

  do i = 1, nNuc

    write(*,*)
    write(*,'(A,I6)') ' Atom index            : ', i
    write(*,'(A,F12.6)') ' Nuclear charge (Z)    : ', ZNuc(i)
    write(*,'(A)') ' Coordinates (x, y, z) : '
    write(*,'(3F18.10)') (rNuc(i,j), j=1,ncart)

  end do

  write(*,*)
  write(*,*) '============================================================'
  write(*,'(A,F20.10,A)') ' Nuclear repulsion energy = ', ENuc,' au'
  write(*,*) '============================================================'
  write(*,*)

end subroutine 
