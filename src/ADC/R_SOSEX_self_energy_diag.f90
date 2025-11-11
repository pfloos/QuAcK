subroutine R_SOSEX_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
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

  integer                       :: i,a,p,m
  double precision              :: num,eps

! Output variables

  double precision,intent(out)  :: Sig(nBas)
  double precision,intent(out)  :: Z(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  Sig(:) = 0d0
  Z(:)   = 0d0

!----------------!
! GW self-energy !
!----------------!
 
! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do m=1,nS

        eps = e(p) - e(i) + Om(m)
        num = 2d0*rhoL(p,i,m)*rhoR(p,i,m)
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = e(p) - e(a) - Om(m)
        num = 2d0*rhoL(p,a,m)*rhoR(p,a,m)
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = e(a) - e(i) + Om(m)
        num = 4d0*rhoL(a,i,m)*rhoR(a,i,m)
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
