double precision function RGW_SRG_ReSigC(p,w,s,nBas,nC,nO,nV,nR,nS,e,Om,rho)

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

! Local variables

  integer                       :: i,a,m
  double precision              :: Dpim,Dpam

! Initialize 

  RGW_SRG_ReSigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do m=1,nS
      Dpim = w - e(i) + Om(m)
      RGW_SRG_ReSigC = RGW_SRG_ReSigC &
                     + 2d0*rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do m=1,nS
      Dpam = w - e(a) - Om(m)
      RGW_SRG_ReSigC = RGW_SRG_ReSigC & 
                     + 2d0*rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam
    end do
  end do

end function 
