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
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs),rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Omega2s(nOOs),Omega2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs),rho2t(nBas,nBas,nOOt)

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
      cd = 0
      do c=nO+1,nBas-nR
        do d=nO+1,c
          cd = cd + 1
          eps = e(p) + e(i) - Omega1s(cd)
          Z(p) = Z(p) - 2d0*rho1s(p,i,cd)**2/eps**2
        enddo
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      kl = 0
      do k=nC+1,nO
        do l=nC+1,k
          kl = kl + 1
          eps = e(p) + e(a) - Omega2s(kl)
          Z(p) = Z(p) - 2d0*rho2s(p,a,kl)**2/eps**2
        enddo
      enddo
    enddo
  enddo

!----------------------------------------------
! Triplet part of the T-matrix self-energy
!----------------------------------------------

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      cd = 0
      do c=nO+1,nBas-nR
        do d=nO+1,c-1
          cd = cd + 1
          eps = e(p) + e(i) - Omega1t(cd)
          Z(p) = Z(p) - 2d0*rho1t(p,i,cd)**2/eps**2
        enddo
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      kl = 0
      do k=nC+1,nO
        do l=nC+1,k-1
          kl = kl + 1
          eps = e(p) + e(a) - Omega2t(kl)
          Z(p) = Z(p) - 2d0*rho2t(p,a,kl)**2/eps**2
        enddo
      enddo
    enddo
  enddo

! Compute renormalization factor from derivative of SigT
 
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine renormalization_factor_Tmatrix
