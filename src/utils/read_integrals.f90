subroutine read_integrals(nBas_AOs, S, T, V, Hc, G)

! Read one- and two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas_AOs

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

! Open file with integrals

  debug = .false.

  lambda = 1d0

  print*, 'Scaling integrals by ',lambda

  open(unit=8 ,file='int/Ov.dat')
  open(unit=9 ,file='int/Kin.dat')
  open(unit=10,file='int/Nuc.dat')

  open(unit=21,file='int/x.dat')
  open(unit=22,file='int/y.dat')
  open(unit=23,file='int/z.dat')

! Read overlap integrals

  S(:,:) = 0d0
  do 
    read(8,*,end=8) mu,nu,Ov
    S(mu,nu) = Ov
    S(nu,mu) = Ov
  end do
  8 close(unit=8)

! Read kinetic integrals

  T(:,:) = 0d0
  do 
    read(9,*,end=9) mu,nu,Kin
    T(mu,nu) = Kin
    T(nu,mu) = Kin
  end do
  9 close(unit=9)

! Read nuclear integrals

  V(:,:) = 0d0
  do 
    read(10,*,end=10) mu,nu,Nuc
    V(mu,nu) = Nuc
    V(nu,mu) = Nuc
  end do
  10 close(unit=10)

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
  open(unit=11, file='int/ERI.bin', form='unformatted', access='stream')
  read(11) G
  close(11)


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
