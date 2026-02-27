double precision function CVS_UGW_Re_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,nCVS,nFC,occupations,virtuals,e,Om,rho)

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
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  integer,intent(in)            :: nCVS,nFC
  integer,intent(in)            :: occupations(nO-nFC)
  integer,intent(in)            :: virtuals(nBas-nO)

! Local variables

  integer                       :: i,a,m
  double precision              :: num,eps

! Initialize 

  CVS_UGW_Re_dSigC = 0d0

! Occupied part of the correlation self-energy

  do i=1,nO-nFC
    do m=1,nS
      eps = w - e(occupations(i)) + Om(m)
      num = rho(p,occupations(i),m)**2
      CVS_UGW_Re_dSigC = CVS_UGW_Re_dSigC - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nCVS+1,nBas-nO
    do m=1,nS
      eps = w - e(virtuals(a)) - Om(m)
      num = rho(p,virtuals(a),m)**2
      CVS_UGW_Re_dSigC = CVS_UGW_Re_dSigC - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do
  end do

end function 
