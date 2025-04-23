subroutine complex_RGW_SRG_SigC_dSigC(flow,p,eta,nBas,nC,nO,nV,nR,nS,Re_w,Im_w,Re_e,Im_e,Om,rho,Re_SigC,Im_SigC,Re_DS,Im_DS)

 ! Complute diagonal of the correlation part of the self-energy and its derivative fully complex

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
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
  double precision              :: eps,s
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
  s = flow

! Compute self energy and its derivative

! Occupied part

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(p,rho,eta,nS,nC,nO,nR,Re_w,Im_w,Re_e,Im_e,Om,s), &
 !$OMP PRIVATE(m,i,eps,num,eta_tilde,tmp) &
 !$OMP REDUCTION(+:Re_SigC,Im_SigC,Re_DS,Im_DS)
  do i=nC+1,nO
    do m=1,nS
      eps = Re_w - Re_e(i) + real(Om(m))
      eta_tilde = eta - Im_w + Im_e(i)  - aimag(Om(m))
      num = 2d0*rho(p,i,m)**2*(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))

      tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),&
              eta_tilde/(eps**2 + eta_tilde**2),kind=8)
      Re_SigC = Re_SigC + real(tmp)
      Im_SigC = Im_SigC + aimag(tmp)

      tmp = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
              -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
      Re_DS    = Re_DS   + real(tmp)
      Im_DS   = Im_DS    + aimag(tmp) 
    end do
  end do
  !$OMP END PARALLEL DO

! Virtual part
  !$OMP PARALLEL &
  !$OMP SHARED(p,nBas,rho,eta,nS,nC,nO,nR,Re_w,Im_w,Re_e,Im_e,Om,s) &
  !$OMP PRIVATE(m,a,eps,tmp,eta_tilde,num) &
  !$OMP REDUCTION(+:Re_SigC,Im_SigC,Re_DS,Im_DS) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do a=nO+1,nBas-nR
    do m=1,nS
      eps = Re_w - Re_e(a) - real(Om(m))
      eta_tilde = eta + Im_w  - Im_e(a) - aimag(Om(m))
      num = 2d0*rho(p,a,m)**2*(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))

      tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),-eta_tilde/(eps**2 + eta_tilde**2),kind=8)
      Re_SigC = Re_SigC + real(tmp)
      Im_SigC = Im_SigC + aimag(tmp)

      tmp = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
              2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
      Re_DS    = Re_DS   + real(tmp)
      Im_DS   = Im_DS    + aimag(tmp) 
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine 
