subroutine regularized_renormalization_factor(COHSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,Z)

! Compute the regularized version of the GW renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
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

  integer                       :: i,a,p,jb
  double precision              :: eps

  double precision              :: kappa
  double precision              :: fk,dfk

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Initialize

  Z(:)  = 0d0

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  kappa = 1d0

! static COHSEX approximation

  if(COHSEX) then
    
    Z(:) = 1d0
    return
 
  else

    ! Occupied part of the correlation self-energy

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do jb=1,nS
          eps = e(p) - e(i) + Omega(jb) 
          fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
          dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)
          Z(p) = Z(p)  - 2d0*rho(p,i,jb)**2*dfk
        end do
      end do
    end do

    ! Virtual part of the correlation self-energy

    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(p) - e(a) - Omega(jb) 
          fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
          dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)
          Z(p) = Z(p)  - 2d0*rho(p,a,jb)**2*dfk
        end do
      end do
    end do

  end if

! Compute renormalization factor from derivative of SigC
 
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine regularized_renormalization_factor
