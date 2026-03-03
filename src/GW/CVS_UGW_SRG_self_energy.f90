subroutine CVS_UGW_SRG_self_energy(flow,nBas,nC,nO,nV,nR,nS,nCVS,nFC,occupations,virtuals,e,Om,rho,EcGM,SigC,Z)

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
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)

! Local variables

  integer                       :: ispin
  integer                       :: i,j,a,b
  integer                       :: p,q
  integer                       :: m
  double precision              :: Dpim,Dqim,Dpam,Dqam,Diam
  double precision              :: s
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)

! SRG flow parameter 

  s = flow

! Initialize 

  SigC(:,:,:) = 0d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nBas,nR,e,Om,nFC,occupations) &
  !$OMP PRIVATE(ispin,m,i,q,p,Dpim,Dqim) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do ispin=1,nspin
    do q=1,nBas
      do p=1,nBas
        do m=1,nS
          do i=1,nO(ispin)-nFC(ispin)

            Dpim = e(p,ispin) - e(occupations(i,ispin),ispin) + Om(m)
            Dqim = e(q,ispin) - e(occupations(i,ispin),ispin) + Om(m)
            SigC(p,q,ispin) = SigC(p,q,ispin) & 
                 + rho(p,occupations(i,ispin),m,ispin)*rho(q,occupations(i,ispin),m,ispin)*(1d0-dexp(-s*Dpim*Dpim)*dexp(-s*Dqim*Dqim)) &
                 *(Dpim + Dqim)/(Dpim*Dpim + Dqim*Dqim)

          end do
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nR,nBas,e,Om,nCVS,virtuals) &
  !$OMP PRIVATE(ispin,m,a,q,p,Dpam,Dqam) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do ispin=1,nspin
    do q=1,nBas
      do p=1,nBas
        do m=1,nS
          do a=nCVS(ispin)+1,nBas-nO(ispin)

            Dpam = e(p,ispin) - e(virtuals(a,ispin),ispin) - Om(m)
            Dqam = e(q,ispin) - e(virtuals(a,ispin),ispin) - Om(m)
            SigC(p,q,ispin) = SigC(p,q,ispin) &
                 + rho(p,virtuals(a,ispin),m,ispin)*rho(q,virtuals(a,ispin),m,ispin)*(1d0-exp(-s*Dpam*Dpam)*exp(-s*Dqam*Dqam)) &
                 *(Dpam + Dqam)/(Dpam*Dpam + Dqam*Dqam)

          end do
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
    do p=1,nBas
      do i=1,nO(ispin)-nFC(ispin)
        do m=1,nS
          Dpim = e(p,ispin) - e(occupations(i,ispin),ispin) + Om(m)
          Z(p,ispin) = Z(p,ispin) - rho(p,occupations(i,ispin),m,ispin)**2*(1d0-dexp(-2d0*s*Dpim*Dpim))/Dpim**2
        end do
      end do
    end do
  end do

  ! Virtual part of the renormalization factor

  do ispin=1,nspin
    do p=1,nBas
      do a=nCVS(ispin)+1,nBas-nO(ispin)
        do m=1,nS
          Dpam = e(p,ispin) - e(virtuals(a,ispin),ispin) - Om(m)
          Z(p,ispin) = Z(p,ispin)  - rho(p,virtuals(a,ispin),m,ispin)**2*(1d0-dexp(-2d0*s*Dpam*Dpam))/Dpam**2
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
  do i=1,nO(ispin)-nFC(ispin)
    do a=nCVS(ispin)+1,nBas-nO(ispin)
      do m=1,nS
        Diam = e(virtuals(a,ispin),ispin) - e(occupations(i,ispin),ispin) + Om(m)
        EcGM = EcGM - rho(virtuals(a,ispin),occupations(i,ispin),m,ispin)*rho(virtuals(a,ispin),occupations(i,ispin),m,ispin)*(1d0-exp(-2d0*s*Diam*Diam))/Diam 
      end do
    end do
  end do
  end do

end subroutine 
