subroutine read_integrals(nBas,S,T,V,Hc,G)

! Read one- and two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,la,si
  double precision              :: Ov,Kin,Nuc,ERI
  double precision              :: lambda

! Output variables

  double precision,intent(out)  :: S(nBas,nBas)
  double precision,intent(out)  :: T(nBas,nBas)
  double precision,intent(out)  :: V(nBas,nBas)
  double precision,intent(out)  :: Hc(nBas,nBas)
  double precision,intent(out)  :: G(nBas,nBas,nBas,nBas)

! Open file with integrals

  debug = .false.

  lambda = 1d0

  print*, 'Scaling integrals by ',lambda

  open(unit=8 ,file='int/Ov.dat')
  open(unit=9 ,file='int/Kin.dat')
  open(unit=10,file='int/Nuc.dat')
  open(unit=11,file='int/ERI.dat')

! Read overlap integrals

  S(:,:) = 0d0
  do 
    read(8,*,end=8) mu,nu,Ov
    S(mu,nu) = Ov
    S(nu,mu) = Ov
  enddo
  8 close(unit=8)

! Read kinetic integrals

  T(:,:) = 0d0
  do 
    read(9,*,end=9) mu,nu,Kin
    T(mu,nu) = Kin
    T(nu,mu) = Kin
  enddo
  9 close(unit=9)

! Read nuclear integrals

  V(:,:) = 0d0
  do 
    read(10,*,end=10) mu,nu,Nuc
    V(mu,nu) = Nuc
    V(nu,mu) = Nuc
  enddo
  10 close(unit=10)

! Define core Hamiltonian

  Hc(:,:) = T(:,:) + V(:,:)

! Read nuclear integrals

  G(:,:,:,:) = 0d0
  do 
    read(11,*,end=11) mu,nu,la,si,ERI

    ERI = lambda*ERI
!   <12|34>
    G(mu,nu,la,si) = ERI
!   <32|14>
    G(la,nu,mu,si) = ERI
!   <14|32>
    G(mu,si,la,nu) = ERI
!   <34|12>
    G(la,si,mu,nu) = ERI
!   <41|23>
    G(si,mu,nu,la) = ERI
!   <23|41>
    G(nu,la,si,mu) = ERI
!   <21|43>
    G(nu,mu,si,la) = ERI
!   <43|21>
    G(si,la,nu,mu) = ERI
  enddo
  11 close(unit=11)


! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Overlap integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,S)
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Kinetic integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,T)
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Nuclear integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,V)
    write(*,*)
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Electron repulsion integrals'
    write(*,'(A28)') '----------------------'
    do la=1,nBas
      do si=1,nBas
        call matout(nBas,nBas,G(1,1,la,si))
      enddo
    enddo
    write(*,*)
  endif

end subroutine read_integrals
