subroutine self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,              & 
                                    Omega1s,rho1s,Omega2s,rho2s,Omega1t,rho1t,Omega2t,rho2t, & 
                                    SigT)

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

  double precision,intent(out)  :: SigT(nBas)

! Initialize 

  SigT(:) = 0d0

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
          SigT(p) = SigT(p) + 2d0*rho1s(p,i,cd)**2/eps
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
          SigT(p) = SigT(p) + 2d0*rho2s(p,a,kl)**2/eps
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
          SigT(p) = SigT(p) + 2d0*rho1t(p,i,cd)**2/eps
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
          SigT(p) = SigT(p) + 2d0*rho2t(p,a,kl)**2/eps
        enddo
      enddo
    enddo
  enddo

end subroutine self_energy_Tmatrix_diag
