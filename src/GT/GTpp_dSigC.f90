double precision function GTpp_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: i,a,cd,kl
  double precision              :: eps

! Initialize 

  GTpp_dSigC = 0d0

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  do i=nC+1,nO
    do cd=1,nVV

      eps = w + e(i) - Om1(cd)
      GTpp_dSigC=  GTpp_dSigC - rho1(p,i,cd)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2

    enddo
  enddo

!----------------------------------------------
! Virtual part of the T-matrix self-energy
!
!----------------------------------------------

  do a=nO+1,nBas-nR
    do kl=1,nOO
      eps = w + e(a) - Om2(kl)
      GTpp_dSigC = GTpp_dSigC - rho2(p,a,kl)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    enddo
  enddo
     
end function 
