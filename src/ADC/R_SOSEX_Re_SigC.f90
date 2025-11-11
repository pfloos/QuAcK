double precision function R_SOSEX_Re_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR)

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
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,m
  double precision              :: num,eps

! Initialize 

  R_SOSEX_Re_SigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do m=1,nS
      eps = w - e(i) + Om(m)
      num = 2d0*rhoL(p,i,m)*rhoR(p,i,m)
      R_SOSEX_Re_SigC = R_SOSEX_Re_SigC + num*eps/(eps**2 + eta**2)
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do m=1,nS
      eps = w - e(a) - Om(m)
      num = 2d0*rhoL(p,a,m)*rhoR(p,a,m)
      R_SOSEX_Re_SigC = R_SOSEX_Re_SigC + num*eps/(eps**2 + eta**2)
    end do
  end do

end function 
