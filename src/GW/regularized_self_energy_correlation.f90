subroutine regularized_self_energy_correlation(COHSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,SigC)

! Compute correlation part of the regularized self-energy

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

  integer                       :: i,j,a,b
  integer                       :: p,q,r
  integer                       :: jb
  double precision              :: eps

  double precision              :: kappa
  double precision              :: fk

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:,:) = 0d0

!---------------------------------------------!
! Parameters for regularized MP2 calculations !
!---------------------------------------------!

  kappa = 1.1d0

!-----------------------------!
! COHSEX static approximation !
!-----------------------------!

  if(COHSEX) then

   ! COHSEX: SEX of the COHSEX correlation self-energy

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do i=nC+1,nO
          do jb=1,nS
            SigC(p,q) = SigC(p,q) + 4d0*rho(p,i,jb)*rho(q,i,jb)/Omega(jb)
          end do
        end do
      end do
    end do
 
    ! COHSEX: COH part of the COHSEX correlation self-energy
 
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do r=nC+1,nBas-nR
          do jb=1,nS
            SigC(p,q) = SigC(p,q) - 2d0*rho(p,r,jb)*rho(q,r,jb)/Omega(jb)
          end do
        end do
      end do
    end do

    EcGM = 0d0
    do i=nC+1,nO
      EcGM = EcGM + 0.5d0*SigC(i,i)
    end do

  else

!----------------!
! GW self-energy !
!----------------!

  ! Occupied part of the correlation self-energy

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do i=nC+1,nO
          do jb=1,nS
            eps = e(p) - e(i) + Omega(jb)
            fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
            SigC(p,q) = SigC(p,q) + 2d0*rho(p,i,jb)*rho(q,i,jb)*fk
          end do
        end do
      end do
    end do
 
    ! Virtual part of the correlation self-energy
 
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do a=nO+1,nBas-nR
          do jb=1,nS
            eps = e(p) - e(a) - Omega(jb)
            fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
            SigC(p,q) = SigC(p,q) + 2d0*rho(p,a,jb)*rho(q,a,jb)*fk
          end do
        end do
      end do
    end do

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(a) - e(i) + Omega(jb)
          fk  = (1d0 - exp(-kappa*abs(eps)))**2/eps
          EcGM = EcGM - 4d0*rho(a,i,jb)**2*fk
        end do
      end do
    end do

  end if

end subroutine regularized_self_energy_correlation
