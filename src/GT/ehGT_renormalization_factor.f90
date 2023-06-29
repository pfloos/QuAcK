subroutine ehGT_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR,Z)

! Compute renormalization factor for GW

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
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

! Local variables

  integer                       :: p,i,a,m
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Initialize

  Z(:)  = 0d0
 
! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do m=1,nS
        eps = e(p) - e(i) + Om(m) 
        Z(p) = Z(p) - rhoL(p,i,m)*rhoR(p,i,m)*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS
        eps = e(p) - e(a) - Om(m) 
        Z(p) = Z(p) - rhoL(p,a,m)*rhoR(p,a,m)*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

! Compute renormalization factor from derivative of SigC
 
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
