subroutine USRG_self_energy(eta,nBas,nC,nO,nV,nR,nS,e,Om,rho,EcGM,SigC,Z)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
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
  integer                       :: p,q,r
  integer                       :: m
  double precision              :: Dpim,Dqim,Dpam,Dqam,Diam
  double precision              :: t1,t2
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)

! Initialize 

  SigC(:,:,:) = 0d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

  ! Occupied part of the correlation self-energy

  call wall_time(t1)

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,eta,nS,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(ispin,m,i,q,p,Dpim,Dqim) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do ispin=1,nspin
  do q=nC(ispin)+1,nBas-nR(ispin)
     do p=nC(ispin)+1,nBas-nR(ispin)
        do m=1,nS
           do i=nC(ispin)+1,nO(ispin)
              Dpim = e(p,ispin) - e(i,ispin) + Om(m)
              Dqim = e(q,ispin) - e(i,ispin) + Om(m)
              SigC(p,q,ispin) = SigC(p,q,ispin) & 
                   + rho(p,i,m,ispin)*rho(q,i,m,ispin)*(1d0-dexp(-eta*Dpim*Dpim)*dexp(-eta*Dqim*Dqim)) &
                   *(Dpim + Dqim)/(Dpim*Dpim + Dqim*Dqim)
           end do
        end do
     end do
  end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

 call wall_time(t2)
 print *, "first loop", (t2-t1)

! Virtual part of the correlation self-energy

 call wall_time(t1)
 !$OMP PARALLEL &
 !$OMP SHARED(SigC,rho,eta,nS,nC,nO,nR,nBas,e,Om) &
 !$OMP PRIVATE(ispin,m,a,q,p,Dpam,Dqam) &
 !$OMP DEFAULT(NONE)
 !$OMP DO
 do ispin=1,nspin
 do q=nC(ispin)+1,nBas-nR(ispin)
    do p=nC(ispin)+1,nBas-nR(ispin)
       do m=1,nS
          do a=nO(ispin)+1,nBas-nR(ispin)
             Dpam = e(p,ispin) - e(a,ispin) - Om(m)
             Dqam = e(q,ispin) - e(a,ispin) - Om(m)
             SigC(p,q,ispin) = SigC(p,q,ispin) &
                  + rho(p,a,m,ispin)*rho(q,a,m,ispin)*(1d0-exp(-eta*Dpam*Dpam)*exp(-eta*Dqam*Dqam)) &
                  *(Dpam + Dqam)/(Dpam*Dpam + Dqam*Dqam)
          end do
       end do
    end do
 end do
 end do
 !$OMP END DO
 !$OMP END PARALLEL

 call wall_time(t2)
  print *, "second loop", (t2-t1)
 

! Initialize

  Z(:,:)  = 0d0

  do ispin=1,nspin
  do p=nC(ispin)+1,nBas-nR(ispin)
     do i=nC(ispin)+1,nO(ispin)
        do m=1,nS
           Dpim = e(p,ispin) - e(i,ispin) + Om(m)
           Z(p,ispin) = Z(p,ispin) - rho(p,i,m,ispin)**2*(1d0-dexp(-2d0*eta*Dpim*Dpim))/Dpim**2
        end do
     end do
  end do
  end do

  ! Virtual part of the correlation self-energy

  do ispin=1,nspin
  do p=nC(ispin)+1,nBas-nR(ispin)
     do a=nO(ispin)+1,nBas-nR(ispin)
        do m=1,nS
           Dpam = e(p,ispin) - e(a,ispin) - Om(m)
           Z(p,ispin) = Z(p,ispin)  - rho(p,a,m,ispin)**2*(1d0-dexp(-2d0*eta*Dpam*Dpam))/Dpam**2
        end do
     end do
  end do
  end do

! Compute renormalization factor from derivative of SigC

  Z(:,:) = 1d0/(1d0 - Z(:,:))

! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do ispin=1,nspin
  do i=nC(ispin)+1,nO(ispin)
    do a=nO(ispin)+1,nBas-nR(ispin)
      do m=1,nS
        Diam = e(a,ispin) - e(i,ispin) + Om(m)
        EcGM = EcGM - rho(a,i,m,ispin)*rho(a,i,m,ispin)*(1d0-exp(-2d0*eta*Diam*Diam))/Diam 
      end do
    end do
  end do
  end do

end subroutine 
