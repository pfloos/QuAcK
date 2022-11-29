subroutine self_energy_correlation_SRG(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,SigC)

! Compute correlation part of the self-energy

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
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,r
  integer                       :: m
  double precision              :: Dpim,Dqim,Dpam,Dqam

! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nBas)

! Initialize 

  SigC(:,:) = 0d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do m=1,nS
          Dpim = e(p) - e(i) + Omega(m)
          Dqim = e(q) - e(i) + Omega(m)
          SigC(p,q) = SigC(p,q) + 2d0*rho(p,i,m)*rho(q,i,m)*(Dpim + Dqim)/(Dpim**2 + Dqim**2)
        end do
      end do
    end do
  end do

! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do m=1,nS
          Dpam = e(p) - e(a) - Omega(m)
          Dqam = e(q) - e(a) - Omega(m)
          SigC(p,q) = SigC(p,q) + 2d0*rho(p,a,m)*rho(q,a,m)*(Dpam + Dqam)/(Dpam**2 + Dqam**2)
        end do
      end do
    end do
  end do

! Galitskii-Migdal correlation energy

  EcGM = 0d0
! do i=nC+1,nO
!   do a=nO+1,nBas-nR
!     do m=1,nS
!       EcGM = EcGM - 4d0*rho(a,i,jb)*rho(a,i,jb)*eps/(eps**2 + eta**2)
!     end do
!   end do
! end do

end subroutine self_energy_correlation_SRG
