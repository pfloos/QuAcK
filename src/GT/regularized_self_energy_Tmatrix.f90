subroutine regularized_self_energy_Tmatrix(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Omega1,rho1,Omega2,rho2,SigT)

! Compute the correlation part of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: i,a,p,q,cd,kl
  double precision              :: eps

  double precision              :: kappa
  double precision              :: fk

! Output variables

  double precision,intent(inout)  :: SigT(nBas,nBas)

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  kappa = 1d0

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do cd=1,nVV

          eps = e(p) + e(i) - Omega1(cd)
          fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps

          SigT(p,q) = SigT(p,q) + rho1(p,i,cd)*rho1(q,i,cd)*fk

        enddo
      enddo
    enddo
  enddo

!----------------------------------------------
  ! Virtual part of the T-matrix self-energy
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do kl=1,nOO

          eps = e(p) + e(a) - Omega2(kl)
          fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps

          SigT(p,q) = SigT(p,q) + rho2(p,a,kl)*rho2(q,a,kl)*fk

        enddo
      enddo
    enddo
  enddo

end subroutine regularized_self_energy_Tmatrix
