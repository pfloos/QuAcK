subroutine self_energy_correlation_SOSEX_diag(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,SigC)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS,nspin)
  double precision,intent(in)   :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: i,a,p,q,jb
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0

!-----------------------------
! SOSEX self-energy
!-----------------------------
 
  ! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do jb=1,nS
        do ispin=1,nspin
          eps = e(p) - e(i) + Omega(jb,ispin)
          SigC(p) = SigC(p) + 2d0*rho(p,i,jb,ispin)**2*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do jb=1,nS
        do ispin=1,nspin
          eps = e(p) - e(a) - Omega(jb,ispin)
          SigC(p) = SigC(p) + 2d0*rho(p,a,jb,ispin)**2*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

  ! GM correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do jb=1,nS
        do ispin=1,nspin
          eps = e(a) - e(i) + Omega(jb,ispin)
          EcGM = EcGM - 4d0*rho(a,i,jb,ispin)*rho(a,i,jb,ispin)*eps/(eps**2 + eta**2)
        end do
      end do
    end do
  end do

end subroutine self_energy_correlation_SOSEX_diag
