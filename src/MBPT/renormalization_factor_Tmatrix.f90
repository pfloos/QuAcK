subroutine renormalization_factor_Tmatrix(alpha,eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Omega1,rho1,Omega2,rho2,Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: alpha
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

  integer                       :: i,j,k,l,a,b,c,d,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do cd=1,nVV
        eps  = e(p) + e(i) - Omega1(cd)
        Z(p) = Z(p) - (rho1(p,i,cd)/eps)**2
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOO
        eps  = e(p) + e(nO+a) - Omega2(kl)
        Z(p) = Z(p) - (rho2(p,nO+a,kl)/eps)**2
      enddo
    enddo
  enddo

end subroutine renormalization_factor_Tmatrix
