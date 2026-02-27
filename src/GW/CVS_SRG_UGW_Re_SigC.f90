double precision function CVS_UGW_SRG_Re_SigC(p,w,s,nBas,nC,nO,nV,nR,nS,nCVS,nFC,occupations,virtuals,e,Om,rho)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: s
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
  double precision              :: Dpim,Dpam

! Initialize 

  CVS_UGW_SRG_Re_SigC = 0d0

! Occupied part of the correlation self-energy

  do i=1,nO-nFC
    do m=1,nS
      Dpim = w - e(occupations(i)) + Om(m)
      CVS_UGW_SRG_Re_SigC = CVS_UGW_SRG_Re_SigC &
                      + rho(p,occupations(i),m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nCVS+1,nBas-nO
    do m=1,nS
      Dpam = w - e(virtuals(a)) - Om(m)
      CVS_UGW_SRG_Re_SigC = CVS_UGW_SRG_Re_SigC & 
                      + rho(p,virtuals(a),m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam
    end do
  end do

end function 
