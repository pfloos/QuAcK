subroutine read_F12_integrals(nBas,S,C,F,Y,FC)

! Read one- and two-electron integrals from files

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  integer                       :: mu,nu,la,si,ka,ta
  double precision              :: ERI,F12,Yuk,F13C12,ExpS

! Output variables

  double precision,intent(out)  :: C(nBas,nBas,nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas,nBas,nBas)
  double precision,intent(out)  :: Y(nBas,nBas,nBas,nBas)
  double precision,intent(out)  :: FC(nBas,nBas,nBas,nBas,nBas,nBas)

  debug = .false.

! Open file with integrals

  open(unit=21,file='int/ERI.dat')
  open(unit=22,file='int/F12.dat')
  open(unit=23,file='int/Yuk.dat')
  open(unit=31,file='int/3eInt_Type1.dat')

! Read 1/r12 integrals

  C = 0d0
  do 
    read(21,*,end=21) mu,nu,la,si,ERI
!   <12|34>
    C(mu,nu,la,si) = ERI
!   <32|14>
    C(la,nu,mu,si) = ERI
!   <14|32>
    C(mu,si,la,nu) = ERI
!   <34|12>
    C(la,si,mu,nu) = ERI
!   <41|23>
    C(si,mu,nu,la) = ERI
!   <23|41>
    C(nu,la,si,mu) = ERI
!   <21|43>
    C(nu,mu,si,la) = ERI
!   <43|21>
    C(si,la,nu,mu) = ERI
  enddo
  21 close(unit=21)

! Read f12 integrals

  F = 0d0
  do 
    read(22,*,end=22) mu,nu,la,si,F12
!   <12|34>
    F(mu,nu,la,si) = F12
!   <32|14>
    F(la,nu,mu,si) = F12
!   <14|32>
    F(mu,si,la,nu) = F12
!   <34|12>
    F(la,si,mu,nu) = F12
!   <41|23>
    F(si,mu,nu,la) = F12
!   <23|41>
    F(nu,la,si,mu) = F12
!   <21|43>
    F(nu,mu,si,la) = F12
!   <43|21>
    F(si,la,nu,mu) = F12
  enddo
  22 close(unit=22)

! Read f12/r12 integrals

  Y = 0d0
  do 
    read(23,*,end=23) mu,nu,la,si,Yuk
!   <12|34>
    Y(mu,nu,la,si) = Yuk
!   <32|14>
    Y(la,nu,mu,si) = Yuk
!   <14|32>
    Y(mu,si,la,nu) = Yuk
!   <34|12>
    Y(la,si,mu,nu) = Yuk
!   <41|23>
    Y(si,mu,nu,la) = Yuk
!   <23|41>
    Y(nu,la,si,mu) = Yuk
!   <21|43>
    Y(nu,mu,si,la) = Yuk
!   <43|21>
    Y(si,la,nu,mu) = Yuk
  enddo
  23 close(unit=23)

! Read f13/r12 integrals

  FC = 0d0
  do 
    read(31,*,end=31) mu,nu,la,si,ka,ta,F13C12
    FC(mu,nu,la,si,ka,ta) = F13C12
  enddo
  31 close(unit=31)

! Print results

  if(debug) then

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Electron repulsion integrals'
    write(*,'(A28)') '----------------------'
    do la=1,nBas
      do si=1,nBas
        call matout(nBas,nBas,C(1,1,la,si))
      enddo
    enddo
    write(*,*)

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'F12 integrals'
    write(*,'(A28)') '----------------------'
    do la=1,nBas
      do si=1,nBas
        call matout(nBas,nBas,F(1,1,la,si))
      enddo
    enddo
    write(*,*)

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Yukawa integrals'
    write(*,'(A28)') '----------------------'
    do la=1,nBas
      do si=1,nBas
        call matout(nBas,nBas,Y(1,1,la,si))
      enddo
    enddo
    write(*,*)

 endif

! Read exponent of Slater geminal
  open(unit=4,file='input/geminal')
  read(4,*) ExpS
  close(unit=4)

! Transform two-electron integrals
  
! do mu=1,nBas
!   do nu=1,nBas
!     do la=1,nBas
!       do si=1,nBas
!         F(mu,nu,la,si) = (S(mu,la)*S(nu,si) - F(mu,nu,la,si))/ExpS
!         Y(mu,nu,la,si) = (C(mu,nu,la,si) - Y(mu,nu,la,si))/ExpS
!       enddo
!     enddo
!   enddo
! enddo

end subroutine read_F12_integrals
