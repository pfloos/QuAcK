subroutine read_LR(nBas,G)

! Read the long-range two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  integer                       :: mu,nu,la,si
  double precision              :: ERI
  double precision              :: lambda

! Output variables

  double precision,intent(out)  :: G(nBas,nBas,nBas,nBas)

! Open file with integrals

  lambda = 1d0

  print*, 'Scaling integrals by ',lambda

  open(unit=11,file='int/ERI_lr.dat')

! Read two-electron integrals

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

end subroutine read_LR
