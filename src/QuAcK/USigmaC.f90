double precision function USigmaC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho)

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
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,jb
  double precision              :: eps

! Initialize 

  USigmaC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do jb=1,nS
      eps = w - e(i) + Omega(jb)
      USigmaC = uSigmaC + rho(p,i,jb)**2*eps/(eps**2 + eta**2)  
    end do
  end do

  do a=nO+1,nBas-nR
    do jb=1,nS
      eps = w - e(a) - Omega(jb)
      USigmaC = USigmaC + rho(p,a,jb)**2*eps/(eps**2 + eta**2)
    end do
  end do

end function USigmaC
