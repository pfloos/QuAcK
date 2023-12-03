subroutine GTeh_regularization(nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR)

! Regularize GTeh excitation densities via SRG

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: e(nBas)
  integer,intent(in)            :: Om(nS)

! Local variables

  integer                       :: p,i,a,m
  double precision              :: s
  double precision              :: kappa
  double precision              :: Dpim,Dpam

! Output variables

  double precision,intent(inout):: rhoL(nBas,nBas,nS)
  double precision,intent(inout):: rhoR(nBas,nBas,nS)

! SRG flow parameter 

  s = 100d0

! Regularize excitation densities

  do p=nC+1,nBas-nR
    do m=1,nS

      do i=nC+1,nO
        Dpim = e(p) - e(i) + Om(m)
        kappa = 1d0 - exp(-Dpim*Dpim*s)
        rhoL(i,p,m) = kappa*rhoL(i,p,m)
        rhoR(i,p,m) = kappa*rhoR(i,p,m)
      end do

      do a=nO+1,nBas-nR
        Dpam = e(p) - e(a) - Om(m)
        kappa = 1d0 - exp(-Dpam*Dpam*s)
        rhoL(p,a,m) = kappa*rhoL(p,a,m)
        rhoR(p,a,m) = kappa*rhoR(p,a,m)
      end do

    end do
  end do

end subroutine 
