subroutine UGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nS,e,Om,rho,EcGM,SigC,Z)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: i,j,a,b
  integer                       :: p
  integer                       :: m
  double precision              :: Dpim,Dpam,Diam
  double precision              :: s
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)

! SRG flow parameter 

  s = flow

! Initialize 

  SigC(:,:) = 0d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(ispin,m,i,p,Dpim) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do ispin=1,nspin
    do p=nC(ispin)+1,nBas-nR(ispin)
      do m=1,nS
        do i=nC(ispin)+1,nO(ispin)

           Dpim = e(p,ispin) - e(i,ispin) + Om(m)
           SigC(p,ispin) = SigC(p,ispin) + rho(p,i,m,ispin)**2*(1d0-dexp(-2d0*s*Dpim*Dpim))/Dpim

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nR,nBas,e,Om) &
  !$OMP PRIVATE(ispin,m,a,p,Dpam) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do ispin=1,nspin
    do p=nC(ispin)+1,nBas-nR(ispin)
      do m=1,nS
        do a=nO(ispin)+1,nBas-nR(ispin)

          Dpam = e(p,ispin) - e(a,ispin) - Om(m)
          SigC(p,ispin) = SigC(p,ispin) + rho(p,a,m,ispin)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam

        end do
      end do
    end do
  end do
 !$OMP END DO
 !$OMP END PARALLEL

!------------------------!
! Renormalization factor !
!------------------------!

  Z(:,:)  = 0d0

  ! Occupied part of the renormalization factor

  do ispin=1,nspin
    do p=nC(ispin)+1,nBas-nR(ispin)
      do i=nC(ispin)+1,nO(ispin)
        do m=1,nS
          Dpim = e(p,ispin) - e(i,ispin) + Om(m)
          Z(p,ispin) = Z(p,ispin) - rho(p,i,m,ispin)**2*(1d0-dexp(-2d0*s*Dpim*Dpim))/Dpim**2
        end do
      end do
    end do
  end do

  ! Virtual part of the renormalization factor

  do ispin=1,nspin
    do p=nC(ispin)+1,nBas-nR(ispin)
      do a=nO(ispin)+1,nBas-nR(ispin)
        do m=1,nS
          Dpam = e(p,ispin) - e(a,ispin) - Om(m)
          Z(p,ispin) = Z(p,ispin)  - rho(p,a,m,ispin)**2*(1d0-dexp(-2d0*s*Dpam*Dpam))/Dpam**2
        end do
      end do
    end do
  end do

  Z(:,:) = 1d0/(1d0 - Z(:,:))

!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  do ispin=1,nspin
  do i=nC(ispin)+1,nO(ispin)
    do a=nO(ispin)+1,nBas-nR(ispin)
      do m=1,nS
        Diam = e(a,ispin) - e(i,ispin) + Om(m)
        EcGM = EcGM - rho(a,i,m,ispin)*rho(a,i,m,ispin)*(1d0-exp(-2d0*s*Diam*Diam))/Diam 
      end do
    end do
  end do
  end do

end subroutine 
