double precision function GTeh_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR)

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
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS,2)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS,2)

! Local variables

  integer                       :: i,a,m
  double precision              :: num,eps

! Initialize 

  GTeh_SigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
     do m=1,nS
        eps = w - e(i) + Om(m)
        num = rhoL(i,p,m,1)*rhoR(i,p,m,2)
        GTeh_SigC = GTeh_SigC + num*eps/(eps**2 + eta**2)
     enddo
  enddo

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
     do m=1,nS
        eps = w - e(a) - Om(m)
        num = rhoL(p,a,m,1)*rhoR(p,a,m,2)
        GTeh_SigC = GTeh_SigC + num*eps/(eps**2 + eta**2)
     enddo
  enddo

end function 
