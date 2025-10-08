subroutine complex_RGW_SRG_self_energy(flow,eta,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,Sig,Z)

! Compute correlation part of the self-energy and the renormalization factor 

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
  complex*16,intent(in)         :: e(nOrb)
  complex*16,intent(in)         :: Om(nS)
  complex*16,intent(in)         :: rho(nOrb,nOrb,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,m
  double precision              :: eps_p,eps_q,eta_tilde_p,eta_tilde_q
  double precision              :: eps,eta_tilde,s
  complex*16                    :: num,tmp
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)
  double precision,allocatable  :: Re_Sig(:,:)
  double precision,allocatable  :: Im_Sig(:,:)
  double precision,allocatable  :: Re_Z(:)
  double precision,allocatable  :: Im_Z(:)

! Output variables

  complex*16,intent(out)        :: EcGM
  complex*16,intent(out)        :: Sig(nOrb,nOrb)
  complex*16,intent(out)        :: Z(nOrb)

!----------------!
! GW self-energy !
!----------------!
  allocate(Re_DS(nOrb),Im_DS(nOrb),Re_Z(nOrb),Im_Z(nOrb),Re_Sig(nOrb,nOrb),Im_Sig(nOrb,nOrb))
  
  Re_Sig(:,:) = 0d0
  Im_Sig(:,:) = 0d0
  Re_DS(:) = 0d0
  Im_DS(:) = 0d0
  s = flow

! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Re_Sig,Im_Sig,rho,eta,nS,nC,nO,nOrb,nR,e,Om,s) &
  !$OMP PRIVATE(m,i,q,p,eps_p,eps_q,num,eta_tilde_p,eta_tilde_q,tmp) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do i=nC+1,nO
          eps_p = real(e(p)) - real(e(i)) + real(Om(m))
          eps_q = real(e(q)) - real(e(i)) + real(Om(m))
          eta_tilde_p = eta  - aimag(e(p)) + aimag(e(i)) - aimag(Om(m)) 
          eta_tilde_q = eta  - aimag(e(q)) + aimag(e(i)) - aimag(Om(m)) 
          num = 2d0*rho(p,i,m)*rho(q,i,m)*(1d0 - exp(-s*(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2)))
          tmp = num*cmplx((eps_p + eps_q)/(eps_p**2 + eps_q**2 + eta_tilde_p**2 + eta_tilde_q**2),&
                  (eta_tilde_p + eta_tilde_q)/(eps_p**2 + eps_q**2 + eta_tilde_p**2 + eta_tilde_q**2),kind=8)
          Re_Sig(p,q) = Re_Sig(p,q) + real(tmp)
          Im_Sig(p,q) = Im_Sig(p,q) + aimag(tmp)
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

 ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Re_Sig,Im_Sig,rho,eta,nS,nC,nO,nOrb,nR,e,Om,s) &
  !$OMP PRIVATE(m,a,q,p,eps_p,eps_q,num,eta_tilde_p,eta_tilde_q,tmp) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do a=nO+1,nOrb-nR
 
          eps_p = real(e(p)) - real(e(a)) - real(Om(m))
          eps_q = real(e(q)) - real(e(a)) - real(Om(m))
          eta_tilde_p = eta  + aimag(e(p)) - aimag(e(a)) - aimag(Om(m))
          eta_tilde_q = eta  + aimag(e(q)) - aimag(e(a)) - aimag(Om(m))
          num = 2d0*rho(p,a,m)*rho(q,a,m)*(1d0 - exp(-s*(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2)))
          tmp = num*cmplx((eps_p +eps_q)/(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2),&
                -(eta_tilde_p + eta_tilde_q)/(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2),kind=8)
          Re_Sig(p,q) = Re_Sig(p,q) + real(tmp)
          Im_Sig(p,q) = Im_Sig(p,q) + aimag(tmp)
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

!------------------------!
! Renormalization factor !
!------------------------!

 !Occupied part of the renormalization factor

 !$OMP PARALLEL &
 !$OMP SHARED(Re_DS,Im_DS,rho,eta,nS,nC,nO,nOrb,nR,e,Om,s) &
 !$OMP PRIVATE(m,i,p,eps,num,eta_tilde,tmp) &
 !$OMP DEFAULT(NONE)
 !$OMP DO
  do p=nC+1,nOrb-nR
     do m=1,nS
        do i=nC+1,nO
          eps = real(e(p)) - real(e(i)) + real(Om(m))
          eta_tilde = eta  - aimag(e(p)) + aimag(e(i)) - aimag(Om(m)) 
          num = 2d0*rho(p,i,m)*rho(p,i,m)*(1d0-exp(-2d0*s*(eps**2+eta_tilde**2)))
          tmp = num*cmplx(-(eps**2-eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                  -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
          Re_DS(p)   = Re_DS(p)   + real(tmp)
          Im_DS(p)   = Im_DS(p)    + aimag(tmp)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the renormalization factor

  !$OMP PARALLEL &
  !$OMP SHARED(Re_DS,Im_DS,rho,eta,nS,nC,nO,nOrb,nR,e,Om,s) &
  !$OMP PRIVATE(m,a,p,eps,num,eta_tilde,tmp) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do p=nC+1,nOrb-nR
    do m=1,nS
      do a=nO+1,nOrb-nR
 
        eps = real(e(p)) - real(e(a)) - real(Om(m))
        eta_tilde = eta  + aimag(e(p)) - aimag(e(a)) - aimag(Om(m))
        num = 2d0*rho(p,a,m)*rho(p,a,m)*(1d0-exp(-2d0*s*(eps**2+eta_tilde**2)))
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
  
  Z = cmplx(Re_Z,Im_Z,kind=8)
  Sig = cmplx(Re_Sig,Im_Sig,kind=8)

  deallocate(Re_DS)
  deallocate(Im_DS)
  deallocate(Re_Z)
  deallocate(Im_Z)
  deallocate(Re_Sig)
  deallocate(Im_Sig)

!!-------------------------------------!
!! Galitskii-Migdal correlation energy !
!!-------------------------------------!
!
  EcGM = 0d0
!  do m=1,nS
!    do a=nO+1,nOrb-nR
!      do i=nC+1,nO
!
!        eps = e(a) - e(i) + Om(m)
!        num = 4d0*rho(a,i,m)*rho(a,i,m)
!        EcGM = EcGM - num*eps/(eps**2 + eta**2)
!
!      end do
!    end do
!  end do
!
end subroutine 
