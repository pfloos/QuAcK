subroutine regularized_renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Omega1,rho1,Omega2,rho2,Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: i,a,p,cd,kl
  double precision              :: eps

  double precision              :: kappa
  double precision              :: fk,dfk

! Output variables

  double precision,intent(out)  :: Z(nBas)

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  kappa = 1d0

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do cd=1,nVV

        eps = e(p) + e(i) - Omega1(cd)

        fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
        dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)
        Z(p) = Z(p) - rho1(p,i,cd)**2*dfk

      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOO

        eps = e(p) + e(nO+a) - Omega2(kl)

        fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
        dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)

        Z(p) = Z(p) - rho2(p,nO+a,kl)**2*dfk

      enddo
    enddo
  enddo

end subroutine regularized_renormalization_factor_Tmatrix
