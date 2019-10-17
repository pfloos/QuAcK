subroutine renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Omega1,rho1,Omega2,rho2,Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
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

! Initialize

  Z(:)  = 0d0

  ! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      cd = 0
      do c=nO+1,nBas-nR
        do d=nO+1,c
          cd = cd + 1
          eps = e(p) + e(i) - Omega1(cd)
          Z(p) = Z(p) - 2d0*rho1(p,i,cd)**2/eps**2
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
          eps = e(p) + e(a) - Omega2(kl)
          Z(p) = Z(p) - 2d0*rho2(p,a,kl)**2/eps**2
        enddo
      enddo
    enddo
  enddo

! Compute renormalization factor from derivative of SigT
 
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine renormalization_factor_Tmatrix
