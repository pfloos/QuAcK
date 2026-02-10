subroutine read_CAP_integrals(working_dir,nBas,W)

! Read one- and two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  character(len=256),intent(in) :: working_dir

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,ios
  double precision              :: wxyz
  character(len=256)            :: file_path

! Output variables

  double precision,intent(out)  :: W(nBas,nBas)

  W(:,:) = 0d0

! Open file with integrals

  debug = .false.
  file_path = trim(working_dir) // '/int/CAP.dat'
  open(unit=31, file=file_path, status='old', iostat=ios)
  if (ios /= 0) then
    return  ! Silently return zero matrix if open fails
  endif

! Read CAP integrals

  do 
    read(31,*,end=31) mu,nu,wxyz
    W(mu,nu) = wxyz
    W(nu,mu) = wxyz
  end do
  31 close(unit=31)

! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'CAP integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,W)
    write(*,*)
  end if

end subroutine 
