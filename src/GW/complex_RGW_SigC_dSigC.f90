subroutine complex_RGW_SigC_dSigC(p,eta,nBas,nC,nO,nV,nR,nS,Re_w,Im_w,Re_e,Im_e,Om,rho,Re_SigC,Im_SigC,Re_DS,Im_DS)

 ! Complute diagonal of the correlation part of the self-energy and its derivative fully complex

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Re_e(nBas)
  double precision,intent(in)   :: Im_e(nBas)
  double precision,intent(in)   :: Re_w
  double precision,intent(in)   :: Im_w
  complex*16,intent(in)         :: Om(nS)
  complex*16,intent(in)         :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,m
  double precision              :: eps
  double precision              :: eta_tilde
  complex*16                    :: num
  complex*16                    :: tmp

! Output variables

  double precision,intent(out)  :: Re_SigC
  double precision,intent(out)  :: Im_SigC
  double precision,intent(out)  :: Re_DS
  double precision,intent(out)  :: Im_DS

! Initialize 
  Re_SigC = 0d0
  Im_SigC = 0d0
  Re_DS    = 0d0
  Im_DS    = 0d0
 

! Compute self energy and its derivative

! Occupied part

  do i=nC+1,nO
    do m=1,nS
      eps = Re_w - Re_e(i) + real(Om(m))
      eta_tilde = eta - Im_w + Im_e(i)  - aimag(Om(m))
      num = 2d0*rho(p,i,m)**2

      tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),eta_tilde/(eps**2 + eta_tilde**2),kind=8)
      Re_SigC = Re_SigC + real(tmp)
      Im_SigC = Im_SigC + aimag(tmp)

      tmp = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
              -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
      Re_DS    = Re_DS   + real(tmp)
      Im_DS   = Im_DS    + aimag(tmp) 
    end do
  end do

! Virtual part

  do a=nO+1,nBas-nR
    do m=1,nS
      eps = Re_w + Re_e(a) - real(Om(m))
      eta_tilde = eta + Im_w  - Im_e(a) - aimag(Om(m))
      num = 2d0*rho(p,a,m)**2

      tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),-eta_tilde/(eps**2 + eta_tilde**2),kind=8)
      Re_SigC = Re_SigC + real(tmp)
      Im_SigC = Im_SigC + aimag(tmp)

      tmp = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
              2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
      Re_DS    = Re_DS   + real(tmp)
      Im_DS   = Im_DS    + aimag(tmp) 
    end do
  end do
end subroutine 
