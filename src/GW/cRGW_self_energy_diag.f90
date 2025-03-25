subroutine cRGW_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,Re_Sig,Im_Sig,Re_Z,Im_Z,CAP)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: CAP(nOrb,nOrb)

! Local variables

  integer                       :: i,a,p,m
  double precision              :: num,eps
  double precision              :: eta_tilde
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)
  
! Output variables

  double precision,intent(out)  :: Re_Sig(nBas)
  double precision,intent(out)  :: Im_Sig(nBas)
  double precision,intent(out)  :: Re_Z(nBas)
  double precision,intent(out)  :: Im_Z(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 
  allocate(Re_DS(nBas),Im_DS(nBas))
  Re_Sig(:) = 0d0
  Im_Sig(:) = 0d0
  Re_DS(:)   = 0d0
  Im_DS(:)   = 0d0

!----------------!
! GW self-energy !
!----------------!
 
! Occupied part of the correlation self-energy
  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do m=1,nS
        eps = e(p) - e(i) + Om(m)
        eta_tilde = eta  - CAP(p,p) + CAP(i,i) 
        num = 2d0*rho(p,i,m)**2
        Re_Sig(p) = Re_Sig(p) + num*eps/(eps**2 + eta_tilde**2)
        Im_Sig(p) = Im_Sig(p) + num*eta_tilde/(eps**2 + eta_tilde**2)
        Re_DS(p)   = Re_DS(p)   - num*(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2
        Im_DS(p)   = Im_DS(p)    - 2*num*eta_tilde*eps/(eps**2 + eta_tilde**2)**2

      end do
    end do
  end do

! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = e(p) - e(a) - Om(m)
        eta_tilde = eta  + CAP(p,p) - CAP(a,a)
        num = 2d0*rho(p,a,m)**2
        Re_Sig(p) = Re_Sig(p) + num*eps/(eps**2 + eta_tilde**2)
        Im_Sig(p) = Im_Sig(p) - num*eta_tilde/(eps**2 + eta_tilde**2)
        Re_DS(p)   = Re_DS(p)   - num*(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2
        Im_DS(p)   = Im_DS(p)    + 2*num*eta_tilde*eps/(eps**2 + eta_tilde**2)**2
      end do
    end do
  end do

! Galitskii-Migdal correlation energy
! MAYBE MODIFY THIS  TERM
  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = e(a) - e(i) + Om(m)
        num = 4d0*rho(a,i,m)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

! Compute renormalization factor from derivative 
  Re_Z(:) = (1d0-Re_DS(:))/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  Im_Z(:) = Im_DS(:)/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  deallocate(Re_DS,Im_DS)
end subroutine 
