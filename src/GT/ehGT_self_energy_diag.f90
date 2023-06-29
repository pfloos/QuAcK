subroutine ehGT_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR,EcGM,SigC)

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
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,p,q,m
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0

!-----------------------------
! GW self-energy
!-----------------------------

  ! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do m=1,nS
        eps = e(p) - e(i) + Om(m)
        SigC(p) = SigC(p) + rhoL(i,p,m)*rhoR(i,p,m)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS
        eps = e(p) - e(a) - Om(m)
        SigC(p) = SigC(p) + rhoL(p,a,m)*rhoR(p,a,m)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  ! GM correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do m=1,nS
        eps = e(a) - e(i) + Om(m)
        EcGM = EcGM - 2d0*rhoL(a,i,m)*rhoR(a,i,m)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

end subroutine 
