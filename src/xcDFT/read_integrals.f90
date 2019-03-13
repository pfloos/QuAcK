subroutine read_integrals(nBas,S,T,V,Hc,G)

! Read one- and two-electron integrals from files

  implicit none

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,la,si
  double precision              :: Ov,Kin,Nuc,ERI

! Output variables

  double precision,intent(out)  :: S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),G(nBas,nBas,nBas,nBas)

! Open file with integrals

  debug = .false.

  open(unit=8 ,file='int/Ov.dat')
  open(unit=9 ,file='int/Kin.dat')
  open(unit=10,file='int/Nuc.dat')
  open(unit=11,file='int/ERI.dat')

! Read overlap integrals

  S = 0d0
  do 
    read(8,*,end=8) mu,nu,Ov
    S(mu,nu) = Ov
  enddo
  8 close(unit=8)

! Read kinetic integrals

  T = 0d0
  do 
    read(9,*,end=9) mu,nu,Kin
    T(mu,nu) = Kin
  enddo
  9 close(unit=9)

! Read nuclear integrals

  V = 0d0
  do 
    read(10,*,end=10) mu,nu,Nuc
    V(mu,nu) = Nuc
  enddo
  10 close(unit=10)

! Define core Hamiltonian

  Hc = T + V

! Read nuclear integrals

  G = 0d0
  do 
    read(11,*,end=11) mu,nu,la,si,ERI
!   (12|34)
    G(mu,nu,la,si) = ERI
!   (21|34)
    G(nu,mu,la,si) = ERI
!   (12|43)
    G(mu,nu,si,la) = ERI
!   (21|43)
    G(nu,mu,si,la) = ERI
!   (34|12)
    G(la,si,mu,nu) = ERI
!   (43|12)
    G(si,la,mu,nu) = ERI
!   (34|21)
    G(la,si,nu,mu) = ERI
!   (43|21)
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
