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
  double precision,intent(in)   :: rho1s(nBas,nO,nVVs),rho1t(nBas,nO,nVVt)
  double precision,intent(in)   :: Omega2s(nOOs),Omega2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nV,nOOs),rho2t(nBas,nV,nOOt)

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
      do cd=1,nVVs
        eps     = e(p)    + e(i) - Omega1s(cd)
        SigT(p) = SigT(p) + rho1s(p,i,cd)**2/eps
      enddo
    enddo
  enddo

  ! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOOs
        eps     = e(p)    + e(nO+a) - Omega2s(kl)
        SigT(p) = SigT(p) + rho2s(p,a,kl)**2/eps
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
        eps     = e(p)    + e(i) - Omega1t(cd)
        SigT(p) = SigT(p) + rho1t(p,i,cd)**2/eps
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC+1,nBas-nR
    do a=1,nV-nR
      do kl=1,nOOt
        eps     = e(p)    + e(nO+a) - Omega2t(kl)
        SigT(p) = SigT(p) + rho2t(p,a,kl)**2/eps
      enddo
    enddo
  enddo

end subroutine self_energy_Tmatrix_diag
