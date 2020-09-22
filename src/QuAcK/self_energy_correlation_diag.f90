subroutine self_energy_correlation_diag(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,rhox,EcGM,SigC)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
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
  double precision,intent(in)   :: rhox(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,p,q,jb
  double precision              :: eps
  double precision,external     :: SigC_dcgw

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0

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
! SOSEX self-energy *BUG*
!-----------------------------

  elseif(SOSEX) then

    ! SOSEX: occupied part of the correlation self-energy

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do jb=1,nS
          eps = e(p) - e(i) + Omega(jb)
          SigC(p) = SigC(p) - rho(p,i,jb)*rhox(p,i,jb)*eps/(eps**2 + eta**2)
        end do
      end do
    end do

    ! SOSEX: virtual part of the correlation self-energy

    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(p) - e(a) - Omega(jb)
          SigC(p) = SigC(p) - rho(p,a,jb)*rhox(p,a,jb)*eps/(eps**2 + eta**2)
        end do
      end do
    end do

    ! GM correlation energy

    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(a) - e(i) + Omega(jb)
          EcGM = EcGM + rho(a,i,jb)*rhox(a,i,jb)*eps/(eps**2 + eta**2)
        end do
      end do
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
          SigC(p) = SigC(p) + 2d0*rho(p,i,jb)**2*eps/(eps**2 + eta**2)
        end do
      end do
    end do
 
    ! Virtual part of the correlation self-energy
 
    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(p) - e(a) - Omega(jb)
          SigC(p) = SigC(p) + 2d0*rho(p,a,jb)**2*eps/(eps**2 + eta**2)
        end do
      end do
    end do

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do jb=1,nS
          eps = e(a) - e(i) + Omega(jb)
          EcGM = EcGM - 4d0*rho(a,i,jb)*rho(a,i,jb)*eps/(eps**2 + eta**2)
        end do
      end do
    end do

  end if

end subroutine self_energy_correlation_diag
