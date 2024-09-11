double precision function UGW_SRG_Re_dSigC(p,w,s,eta,nBas,nC,nO,nV,nR,nS,e,Om,rho)

! Compute the derivative of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: s
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

! Local variables

  integer                       :: i,a,m
  double precision              :: Dpim,Dpam

! Initialize 

  UGW_SRG_Re_dSigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do m=1,nS
      Dpim = w - e(i) + Om(m)
      UGW_SRG_Re_dSigC = UGW_SRG_Re_dSigC &
                       - rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim**2
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do m=1,nS
      Dpam = w - e(a) - Om(m)
      UGW_SRG_Re_dSigC = UGW_SRG_Re_dSigC & 
                       - rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam**2
    end do
  end do

end function 
