subroutine read_2e_integrals(working_dir,nBas_AOs,G)

! Read two-electron integrals from files

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs
  character(len=256),intent(in) :: working_dir

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,la,si
  double precision              :: ERI
  double precision              :: lambda

! Output variables

  double precision,intent(out)  :: G(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)

  integer                       :: ios
  character(len=256)            :: file_path

! Open file with integrals

  debug = .false.

  lambda = 1d0

  print*, 'Scaling integrals by ',lambda


! Read 2e-integrals

!  ! formatted file
!  open(unit=11, file='int/ERI.dat')
!  G(:,:,:,:) = 0d0
!  do 
!    read(11,*,end=11) mu, nu, la, si, ERI
!    ERI = lambda*ERI
!    G(mu,nu,la,si) = ERI    !   <12|34>
!    G(la,nu,mu,si) = ERI    !   <32|14>
!    G(mu,si,la,nu) = ERI    !   <14|32>
!    G(la,si,mu,nu) = ERI    !   <34|12>
!    G(si,mu,nu,la) = ERI    !   <41|23>
!    G(nu,la,si,mu) = ERI    !   <23|41>
!    G(nu,mu,si,la) = ERI    !   <21|43>
!    G(si,la,nu,mu) = ERI    !   <43|21>
!  end do
!  11 close(unit=11)

  ! binary file
  file_path = trim(working_dir) // '/int/ERI.bin'
  open(unit=11, file=file_path, status='old', action='read', form='unformatted', access='stream', iostat=ios)
    if(ios /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      read(11) G
    endif
  close(unit=11)
  G = G * lambda



! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Electron repulsion integrals'
    write(*,'(A28)') '----------------------'
    do la=1,nBas_AOs
      do si=1,nBas_AOs
        call matout(nBas_AOs, nBas_AOs, G(1,1,la,si))
      end do
    end do
    write(*,*)
  end if

end subroutine 

! ---

subroutine read_2e_integrals_hpc(working_dir, ERI_size, ERI_chem)

  implicit none

  character(len=256), intent(in)  :: working_dir
  integer*8,          intent(in)  :: ERI_size
  double precision,   intent(out) :: ERI_chem(ERI_size)

  integer                         :: ios
  character(len=256)              :: file_path

  file_path = trim(working_dir) // '/int/ERI_chem.bin'
  open(unit=11, file=file_path, status='old', action='read', form='unformatted', access='stream', iostat=ios)
    if(ios /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      read(11) ERI_chem
    endif
  close(unit=11)

  return
end subroutine 

! ---

