subroutine self_energy_Tmatrix_diag_so(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Omega1,rho1,Omega2,rho2,SigT)

! Compute diagonal of the correlation part of the T-matrix self-energy

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

  integer                       :: i,j,k,l,a,b,c,d,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: SigT(nBas)

! Initialize 

  SigT(:) = 0d0

!----------------------------------------------
! T-matrix self-energy in the spinorbital basis
!----------------------------------------------

! Occupied part of the T-matrix self-energy 

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1
          eps = e(p) + e(i) - Omega1(cd)
          SigT(p) = SigT(p) + rho1(p,i,cd)**2/eps
        enddo
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1
          eps = e(p) + e(a) - Omega2(kl)
          SigT(p) = SigT(p) + rho2(p,a,kl)**2/eps
        enddo
      enddo
    enddo
  enddo

end subroutine self_energy_Tmatrix_diag_so
