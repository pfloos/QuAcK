subroutine complex_RGW_SRG_self_energy_diag(flow,eta,nBas,nOrb,nC,nO,nV,nR,nS,Re_e,Im_e,Om,rho,EcGM,Re_Sig,Im_Sig,Re_Z,Im_Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Re_e(nBas)
  double precision,intent(in)   :: Im_e(nBas)
  complex*16,intent(in)         :: Om(nS)
  complex*16,intent(in)         :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,p,m
  double precision              :: eps,s
  complex*16                    :: num
  double precision              :: eta_tilde
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)
  complex*16                    :: tmp
  
! Output variables

  double precision,intent(out)  :: Re_Sig(nBas)
  double precision,intent(out)  :: Im_Sig(nBas)
  double precision,intent(out)  :: Re_Z(nBas)
  double precision,intent(out)  :: Im_Z(nBas)
  complex*16,intent(out)        :: EcGM

! Initialize 
  allocate(Re_DS(nBas),Im_DS(nBas))
  Re_Sig(:) = 0d0
  Im_Sig(:) = 0d0
  Re_DS(:)   = 0d0
  Im_DS(:)   = 0d0

  s = flow

!----------------!
! GW self-energy !
!----------------!
 
! Occupied part of the correlation self-energy
 !$OMP PARALLEL &
 !$OMP SHARED(nBas,Re_Sig,Im_Sig,Re_Z,Im_Z,rho,eta,nS,nC,nO,nOrb,nR,Re_e,Im_e,Om,Re_DS,s,Im_DS), &
 !$OMP PRIVATE(m,i,p,eps,num,eta_tilde,tmp) &
 !$OMP DEFAULT(NONE)
 !$OMP DO
  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do m=1,nS
        eps = Re_e(p) - Re_e(i) + real(Om(m))
        eta_tilde = eta  - Im_e(p) + Im_e(i) - aimag(Om(m)) 
        num = 2d0*rho(p,i,m)**2*(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))
        tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),&
                eta_tilde/(eps**2+eta_tilde**2),kind=8)
        Re_Sig(p) = Re_Sig(p) + real(tmp)
        Im_Sig(p) = Im_Sig(p) + aimag(tmp)
        tmp = num*cmplx(-(eps**2-eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
        Re_DS(p)   = Re_DS(p)   + real(tmp)
        Im_DS(p)   = Im_DS(p)    + aimag(tmp)

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the correlation self-energy
  !$OMP PARALLEL &
  !$OMP SHARED(nBas,Re_Sig,Im_Sig,Re_Z,Im_Z,Re_DS,Im_DS,rho,eta,nS,nC,nO,nOrb,nR,Re_e,Im_e,Om,s) &
  !$OMP PRIVATE(m,a,p,eps,tmp,eta_tilde,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = Re_e(p) - Re_e(a) - real(Om(m))
        eta_tilde = eta  + Im_e(p) - Im_e(a) - aimag(Om(m))
        num = 2d0*rho(p,a,m)**2*(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))
        tmp = num*cmplx(eps/(eps**2 + eta_tilde**2),&
                -eta_tilde/(eps**2 + eta_tilde**2),kind=8)
        Re_Sig(p) = Re_Sig(p) + real(tmp)
        Im_Sig(p) = Im_Sig(p) + aimag(tmp)
        tmp = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
        Re_DS(p)   = Re_DS(p)    + real(tmp)
        Im_DS(p)   = Im_DS(p)    + aimag(tmp)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Compute renormalization factor from derivative 
  Re_Z(:) = (1d0-Re_DS(:))/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  Im_Z(:) = Im_DS(:)/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  deallocate(Re_DS,Im_DS)
end subroutine 
