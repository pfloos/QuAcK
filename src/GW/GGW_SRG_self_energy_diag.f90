subroutine GGW_SRG_self_energy_diag(flow,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,SigC,Z)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p
  integer                       :: m
  double precision              :: Dpim,Dpam,Diam
  double precision              :: s
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

!--------------------!
! SRG flow parameter !
!--------------------!

  s = flow

!--------------------!
! SRG-GW self-energy !
!--------------------!

  SigC(:) = 0d0

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,i,p,Dpim) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do p=nC+1,nOrb-nR
    do m=1,nS
      do i=nC+1,nO

        Dpim = e(p) - e(i) + Om(m)
        SigC(p) = SigC(p) + rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nR,nOrb,e,Om) &
  !$OMP PRIVATE(m,a,p,Dpam) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do p=nC+1,nOrb-nR
    do m=1,nS
      do a=nO+1,nOrb-nR
 
         Dpam = e(p) - e(a) - Om(m)
         SigC(p) = SigC(p) + rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam
 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
 
!------------------------!
! Renormalization factor !
!------------------------!

  Z(:)  = 0d0

  ! Occupied part of the renormlization factor

  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do m=1,nS

        Dpim = e(p) - e(i) + Om(m)
        Z(p) = Z(p) - rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim**2

      end do
    end do
  end do

  ! Virtual part of the renormlization factor

  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do m=1,nS
           Dpam = e(p) - e(a) - Om(m)
           Z(p) = Z(p)  - rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam**2
        end do
     end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do m=1,nS

        Diam = e(a) - e(i) + Om(m)
        EcGM = EcGM - rho(a,i,m)*rho(a,i,m)*(1d0-exp(-2d0*s*Diam*Diam))/Diam 

      end do
    end do
  end do

end subroutine 
