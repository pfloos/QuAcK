double precision function dUSigmaC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho)

! Compute the derivative of the correlation part of the self-energy

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
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: i,a,jb
  double precision              :: eps

! Initialize 

  dUSigmaC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do jb=1,nS
      eps = w - e(i) + Omega(jb)
      dUSigmaC = dUSigmaC + rho(p,i,jb,1)**2*(eps/(eps**2 + eta**2))**2
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do jb=1,nS
      eps = w - e(a) - Omega(jb)
      dUSigmaC = dUSigmaC + rho(p,a,jb,1)**2*(eps/(eps**2 + eta**2))**2
    end do
  end do

end function dUSigmaC
