subroutine RGW_SRG_self_energy(nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,SigC,Z)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
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
  integer                       :: p,q
  integer                       :: m
  double precision              :: Dpim,Dqim,Dpam,Dqam,Diam
  double precision              :: renorm
  double precision              :: s
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb,nOrb)
  double precision,intent(out)  :: Z(nOrb)

!--------------------!
! SRG flow parameter !
!--------------------!

  s = 500d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

  SigC(:,:) = 0d0

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,Dpim,Dqim,renorm) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do i=nC+1,nO

          Dpim = e(p) - e(i) + Om(m)
          Dqim = e(q) - e(i) + Om(m)
          renorm = (1d0-exp(-s*Dpim*Dpim)*exp(-s*Dqim*Dqim))*(Dpim + Dqim)/(Dpim*Dpim + Dqim*Dqim)
          SigC(p,q) = SigC(p,q) + 2d0*rho(p,i,m)*rho(q,i,m)*renorm

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nR,nOrb,e,Om) &
  !$OMP PRIVATE(m,a,q,p,Dpam,Dqam,renorm) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do a=nO+1,nOrb-nR
 
           Dpam = e(p) - e(a) - Om(m)
           Dqam = e(q) - e(a) - Om(m)
           renorm = (1d0-exp(-s*Dpam*Dpam)*exp(-s*Dqam*Dqam))*(Dpam + Dqam)/(Dpam*Dpam + Dqam*Dqam)
           SigC(p,q) = SigC(p,q) + 2d0*rho(p,a,m)*rho(q,a,m)*renorm
 
        end do
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
        Z(p) = Z(p) - 2d0*rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim**2

      end do
    end do
  end do

  ! Virtual part of the renormlization factor

  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do m=1,nS
           Dpam = e(p) - e(a) - Om(m)
           Z(p) = Z(p)  - 2d0*rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam**2
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
        EcGM = EcGM - 4d0*rho(a,i,m)*rho(a,i,m)*(1d0-exp(-2d0*s*Diam*Diam))/Diam 

      end do
    end do
  end do

end subroutine 
