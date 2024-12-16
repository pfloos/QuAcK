subroutine read_1e_integrals(working_dir,nBas_AOs,S,T,V,Hc)

! Read one-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas_AOs
  character(len=256),intent(in) :: working_dir

! Local variables

  logical                       :: debug
  integer                       :: mu,nu
  double precision              :: Ov,Kin,Nuc

! Output variables

  double precision,intent(out)  :: S(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: T(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: V(nBas_AOs,nBas_AOs)
  double precision,intent(out)  :: Hc(nBas_AOs,nBas_AOs)

  integer                       :: ios
  character(len=256)            :: file_path

! Open file with integrals

  debug = .false.


  ! ---

  ! Read overlap integrals
  file_path = trim(working_dir) // '/int/Ov.dat'
  open(unit=8, file=file_path, status='old', action='read', iostat=ios)
    if(ios /= 0) then
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
  open(unit=9, file=file_path, status='old', action='read', iostat=ios)
    if(ios /= 0) then
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
  open(unit=10, file=file_path, status='old', action='read', iostat=ios)
    if(ios /= 0) then
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
  end if

end subroutine 
