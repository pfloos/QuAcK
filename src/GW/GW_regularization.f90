subroutine GW_regularization(nBas,nC,nO,nV,nR,nS,e,Om,rho)

! Regularize GW excitation densities via SRG

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

  double precision,intent(inout):: rho(nBas,nBas,nS)

! SRG flow parameter 

  s = 100d0

! Regularize excitation densities

  do p=nC+1,nBas-nR
    do m=1,nS

      do i=nC+1,nO
        Dpim = e(p) - e(i) + Om(m)
        kappa = 1d0 - exp(-Dpim*Dpim*s)
        rho(p,i,m) = kappa*rho(p,i,m)
      end do

      do a=nO+1,nBas-nR
        Dpam = e(p) - e(a) - Om(m)
        kappa = 1d0 - exp(-Dpam*Dpam*s)
        rho(p,a,m) = kappa*rho(p,a,m)
      end do

    end do
  end do

end subroutine 
