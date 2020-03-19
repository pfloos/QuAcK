subroutine renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,              & 
                                          Omega1s,rho1s,Omega2s,rho2s,Omega1t,rho1t,Omega2t,rho2t, & 
                                          Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  integer,intent(in)            :: nOOs,nOOt
  integer,intent(in)            :: nVVs,nVVt
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega1s(nVVs),Omega1t(nVVt)
  double precision,intent(in)   :: rho1s(nBas,nO,nVVs),rho1t(nBas,nO,nVVt)
  double precision,intent(in)   :: Omega2s(nOOs),Omega2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nV,nOOs),rho2t(nBas,nV,nOOt)

! Local variables

  integer                       :: i,j,k,l,a,b,c,d,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Initialize

  Z(:)  = 0d0

!----------------------------------------------
! Singlet part of the T-matrix self-energy
!----------------------------------------------

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do cd=1,nVVs
        eps = e(p) + e(i) - Omega1s(cd)
        Z(p) = Z(p) + 1d0*(rho1s(p,i,cd)/eps)**2
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOOs
        eps = e(p) + e(nO+a) - Omega2s(kl)
        Z(p) = Z(p) + 1d0*(rho2s(p,a,kl)/eps)**2
      enddo
    enddo
  enddo

!----------------------------------------------
! Triplet part of the T-matrix self-energy
!----------------------------------------------

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do cd=1,nVVt
        eps = e(p) + e(i) - Omega1t(cd)
        Z(p) = Z(p) + 1d0*(rho1t(p,i,cd)/eps)**2
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOOt
        eps = e(p) + e(nO+a) - Omega2t(kl)
        Z(p) = Z(p) + 1d0*(rho2t(p,a,kl)/eps)**2
      enddo
    enddo
  enddo

! Compute renormalization factor from derivative of SigT
 
  Z(:) = 1d0/(1d0 + Z(:))

end subroutine renormalization_factor_Tmatrix
