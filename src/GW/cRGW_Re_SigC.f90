double precision function cRGW_Re_SigC(p,Re_w,Im_w,eta,nBas,nC,nO,nV,nR,nS,Re_e,Im_e,Om,rho)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: Re_w
  double precision,intent(in)   :: Im_w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Re_e(nBas)
  double precision,intent(in)   :: Im_e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,m
  double precision              :: num,eps
  double precision              :: eta_tilde

! Initialize 

  cRGW_Re_SigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do m=1,nS
      eps = Re_w - Re_e(i) + Om(m)
      eta_tilde = eta - Im_w + Im_e(i)
      num = 2d0*rho(p,i,m)**2
      cRGW_Re_SigC = cRGW_Re_SigC + num*eps/(eps**2 + eta_tilde**2)
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do m=1,nS
      eps = Re_w - Re_e(a) - Om(m)
      eta_tilde = eta + Im_w - Im_e(a)
      num = 2d0*rho(p,a,m)**2
      cRGW_Re_SigC = cRGW_Re_SigC + num*eps/(eps**2 + eta_tilde**2)
    end do
  end do

end function 
