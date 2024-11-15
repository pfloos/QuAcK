subroutine read_integrals(working_dir,nBas_AOs,S,T,V,Hc,G)

! Read one- and two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas_AOs
  character(len=256),intent(in) :: working_dir

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,la,si
  double precision              :: Ov,Kin,Nuc,ERI
  double precision              :: lambda

! Output variables

  double precision,intent(out)  :: S(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: T(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: V(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: G(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)

  integer                       :: status, ios
  character(len=256)            :: file_path

! Open file with integrals

  debug = .false.

  lambda = 1d0

  print*, 'Scaling integrals by ',lambda


  ! ---

  ! Read overlap integrals
  file_path = trim(working_dir) // '/int/Ov.dat'
  open(unit=8, file=file_path, status='old', action='read', iostat=status)
    if(status /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      S(:,:) = 0d0
      do 
        read(8,*,iostat=ios) mu,nu,Ov
        if(ios /= 0) exit
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      end do
    endif
  close(unit=8)

  ! ---

  ! Read kinetic integrals
  file_path = trim(working_dir) // '/int/Kin.dat'
  open(unit=9, file=file_path, status='old', action='read', iostat=status)
    if(status /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      T(:,:) = 0d0
      do 
        read(9,*,iostat=ios) mu,nu,Kin
        if(ios /= 0) exit
        T(mu,nu) = Kin
        T(nu,mu) = Kin
      end do
    endif
  close(unit=9)

  ! ---

  ! Read nuclear integrals
  file_path = trim(working_dir) // '/int/Nuc.dat'
  open(unit=10, file=file_path, status='old', action='read', iostat=status)
    if(status /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      V(:,:) = 0d0
      do 
        read(10,*,iostat=ios) mu,nu,Nuc
        if(ios /= 0) exit
        V(mu,nu) = Nuc
        V(nu,mu) = Nuc
      end do
    endif
  close(unit=10)

  ! ---

  ! Define core Hamiltonian
  Hc(:,:) = T(:,:) + V(:,:)

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
  open(unit=11, file=file_path, status='old', action='read', form='unformatted', access='stream', iostat=status)
    if(status /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      read(11) G
    endif
  close(unit=11)



! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Overlap integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas_AOs,nBas_AOs,S)
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Kinetic integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas_AOs,nBas_AOs,T)
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Nuclear integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas_AOs,nBas_AOs,V)
    write(*,*)
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
