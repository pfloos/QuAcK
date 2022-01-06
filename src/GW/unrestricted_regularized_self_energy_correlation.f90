subroutine unrestricted_regularized_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nSt,e,Omega,rho,SigC,EcGM)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: rho(nBas,nBas,nSt,nspin)

! Local variables

  integer                       :: i,a,p,q,jb
  double precision              :: eps

  double precision              :: kappa
  double precision              :: fk

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas,nspin)
  double precision              :: EcGM(nspin)

! Initialize 

  SigC(:,:,:) = 0d0
  EcGM(:)   = 0d0

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  kappa = 1.1d0

!--------------!
! Spin-up part !
!--------------!

  ! Occupied part of the correlation self-energy

  do p=nC(1)+1,nBas-nR(1)
    do q=nC(1)+1,nBas-nR(1)
      do i=nC(1)+1,nO(1)
        do jb=1,nSt
          eps = e(p,1) - e(i,1) + Omega(jb)
          fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
          SigC(p,q,1) = SigC(p,q,1) + rho(p,i,jb,1)*rho(q,i,jb,1)*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC(1)+1,nBas-nR(1)
    do q=nC(1)+1,nBas-nR(1)
      do a=nO(1)+1,nBas-nR(1)
        do jb=1,nSt
          eps = e(p,1) - e(a,1) - Omega(jb)
          fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
          SigC(p,q,1) = SigC(p,q,1) + rho(p,a,jb,1)*rho(q,a,jb,1)*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! GM correlation energy

  do i=nC(1)+1,nO(1)
    do a=nO(1)+1,nBas-nR(1)
      do jb=1,nSt
        eps = e(a,1) - e(i,1) + Omega(jb)
        fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
        EcGM(1) = EcGM(1) - rho(a,i,jb,1)**2*eps/(eps**2 + eta**2)
      end do
    end do
  end do

!----------------!
! Spin-down part !
!----------------!

  ! Occupied part of the correlation self-energy

  do p=nC(2)+1,nBas-nR(2)
    do q=nC(2)+1,nBas-nR(2)
      do i=nC(2)+1,nO(2)
        do jb=1,nSt
          eps = e(p,2) - e(i,2) + Omega(jb)
          fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
          SigC(p,q,2) = SigC(p,q,2) + rho(p,i,jb,2)*rho(q,i,jb,2)*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC(2)+1,nBas-nR(2)
    do q=nC(2)+1,nBas-nR(2)
      do a=nO(2)+1,nBas-nR(2)
        do jb=1,nSt
          eps = e(p,2) - e(a,2) - Omega(jb)
          fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
          SigC(p,q,2) = SigC(p,q,2) + rho(p,a,jb,2)*rho(q,a,jb,2)*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! GM correlation energy

  do i=nC(2)+1,nO(2)
    do a=nO(2)+1,nBas-nR(2)
      do jb=1,nSt
        eps = e(a,2) - e(i,2) + Omega(jb)
        fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
        EcGM(2) = EcGM(2) - rho(a,i,jb,2)**2*eps/(eps**2 + eta**2)
      end do
    end do
  end do

end subroutine unrestricted_regularized_self_energy_correlation
