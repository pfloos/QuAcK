subroutine regularized_self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,SigC)

! Compute diagonal of the correlation part of the regularized self-energy

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,p,q,jb
  double precision              :: eps
  double precision,external     :: SigC_dcgw

  double precision              :: kappa
  double precision              :: fk

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0

!---------------------------------------------!
! Parameters for regularized MP2 calculations !
!---------------------------------------------!

  kappa = 1.1d0

!-----------------------------
! COHSEX static self-energy
!-----------------------------

  if(COHSEX) then

    ! COHSEX: SEX part of the COHSEX correlation self-energy

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do jb=1,nS
          SigC(p) = SigC(p) + 4d0*rho(p,i,jb)**2/Omega(jb)
        end do
      end do
    end do
 
    ! COHSEX: COH part of the COHSEX correlation self-energy
 
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do jb=1,nS
          SigC(p) = SigC(p) - 2d0*rho(p,q,jb)**2/Omega(jb)
        end do
      end do
    end do

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      EcGM = EcGM - SigC(i)
    end do

!-----------------------------
! GW self-energy
!-----------------------------

  else
 
    ! Occupied part of the correlation self-energy
 
    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do jb=1,nS
          eps = e(p) - e(i) + Omega(jb)
          fk  = (1d0 - exp(-kappa*eps))**2/eps
          SigC(p) = SigC(p) + 2d0*rho(p,i,jb)**2*fk
        end do
      end do
    end do
 
    ! Virtual part of the correlation self-energy
 
    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(p) - e(a) - Omega(jb)
          fk  = (1d0 - exp(-kappa*eps))**2/eps
          SigC(p) = SigC(p) + 2d0*rho(p,a,jb)**2*fk
        end do
      end do
    end do

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(a) - e(i) + Omega(jb)
          fk  = (1d0 - exp(-kappa*eps))**2/eps
          EcGM = EcGM - 4d0*rho(a,i,jb)**2*fk
        end do
      end do
    end do

  end if

end subroutine regularized_self_energy_correlation_diag
