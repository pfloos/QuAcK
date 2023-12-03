subroutine read_dipole_integrals(nBas,R)

! Read one-electron integrals related to dipole moment from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: debug = .false.
  integer                       :: mu,nu
  double precision              :: Dip

! Output variables

  double precision,intent(out)  :: R(nBas,nBas,ncart)

! Open file with integrals

  open(unit=21,file='int/x.dat')
  open(unit=22,file='int/y.dat')
  open(unit=23,file='int/z.dat')

! Read (x,y,z) integrals

  R(:,:,:) = 0d0

  do
    read(21,*,end=21) mu,nu,Dip
    R(mu,nu,1) = Dip
    R(nu,mu,1) = Dip
  end do
  21 close(unit=21)

  do
    read(22,*,end=22) mu,nu,Dip
    R(mu,nu,2) = Dip
    R(nu,mu,2) = Dip
  end do
  22 close(unit=22)

  do
    read(23,*,end=23) mu,nu,Dip
    R(mu,nu,3) = Dip
    R(nu,mu,3) = Dip
  end do
  23 close(unit=23)

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

end subroutine read_dipole_integrals
