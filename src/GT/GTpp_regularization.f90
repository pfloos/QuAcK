subroutine GTpp_regularization(nBas,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2)

! Regularize GTpp excitation densities via SRG

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  integer,intent(in)            :: e(nBas)

  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: Om2(nOO)

! Local variables

  integer                       :: p,i,a,m
  double precision              :: s
  double precision              :: kappa
  double precision              :: Dpim,Dpam

! Output variables

  double precision,intent(inout):: rho1(nBas,nBas,nVV)
  double precision,intent(inout):: rho2(nBas,nBas,nOO)

! SRG flow parameter 

  s = 100d0

! Regularize excitation densities

  do p=nC+1,nBas-nR

    do m=1,nVV
      do i=nC+1,nO
        Dpim = e(p) + e(i) - Om1(m)
        kappa  = exp(-Dpim*Dpim*s)
        rho1(p,i,m) = kappa*rho1(p,i,m)
      end do
    end do

    do m=1,nOO
      do a=nO+1,nBas-nR
        Dpam = e(p) + e(a) - Om2(m)
        kappa  = exp(-Dpam*Dpam*s)
        rho2(p,a,m) = kappa*rho2(p,a,m)
      end do
    end do

  enddo

end subroutine 
