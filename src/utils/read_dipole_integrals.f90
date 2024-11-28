subroutine read_dipole_integrals(working_dir,nBas,R)

! Read one-electron integrals related to dipole moment from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  character(len=256),intent(in) :: working_dir

! Local variables

  logical                       :: debug = .false.
  integer                       :: mu,nu
  double precision              :: Dip

! Output variables

  double precision,intent(out)  :: R(nBas,nBas,ncart)

  integer                       :: status, ios
  character(len=256)            :: file_path


! Open file with integrals

  R(:,:,:) = 0d0

  file_path = trim(working_dir) // '/int/x.dat'
  open(unit=21, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop
  
    else

      do
        read(21, '(I5, I5, E25.17)', iostat=ios) mu, nu, Dip
        if(ios /= 0) exit
        R(mu,nu,1) = Dip
        R(nu,mu,1) = Dip
      end do

    endif

  close(unit=21)

  ! ---

  file_path = trim(working_dir) // '/int/y.dat'
  open(unit=22, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop
  
    else

      do
        read(22, '(I5, I5, E25.17)', iostat=ios) mu, nu, Dip
        if(ios /= 0) exit
        R(mu,nu,2) = Dip
        R(nu,mu,2) = Dip
      end do

    endif

  close(unit=22)

  ! ---

  file_path = trim(working_dir) // '/int/z.dat'
  open(unit=23, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop
  
    else

      do
        read(23, '(I5, I5, E25.17)', iostat=ios) mu, nu, Dip
        if(ios /= 0) exit
        R(mu,nu,3) = Dip
        R(nu,mu,3) = Dip
      end do

    endif

  close(unit=23)

  ! ---


! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') '(mu|x|nu) integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,R(:,:,1))
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') '(mu|y|nu) integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,R(:,:,2))
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') '(mu|z|nu) integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,R(:,:,3))
  end if

end subroutine 
