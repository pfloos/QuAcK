subroutine renormalization_factor_SRG(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,Z)

! Compute renormalization factor for GW

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: p,i,a,m
  double precision              :: Dpim,Dpam

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Initialize

  Z(:)  = 0d0
  
  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do m=1,nS
           Dpim = e(p) - e(i) + Omega(m)
           Z(p) = Z(p)  - 2d0*rho(p,i,m)**2*(1d0-dexp(-2d0*eta*Dpim*Dpim))/Dpim**2
        end do
     end do
  end do

  ! Virtual part of the correlation self-energy
  
  do p=nC+1,nBas-nR
     do a=nO+1,nBas-nR
        do m=1,nS
           Dpam = e(p) - e(a) - Omega(m)
           Z(p) = Z(p)  - 2d0*rho(p,a,m)**2*(1d0-dexp(-2d0*eta*Dpam*Dpam))/Dpam**2
        end do
     end do
  end do



! Compute renormalization factor from derivative of SigC
 
  Z(:) = 1d0/(1d0 - Z(:))
  
end subroutine renormalization_factor_SRG
