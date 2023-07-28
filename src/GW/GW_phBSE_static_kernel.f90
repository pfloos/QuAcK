subroutine GW_phBSE_static_kernel(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,W)

! Compute the second-order static BSE kernel for the resonant block (only for singlets!)

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
  double precision,intent(in)   :: lambda

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Om(nS)

  double precision,intent(in)   :: rho(nBas,nBas,nS)


! Local variables

  double precision              :: chi
  integer                       :: p,q,r,s
  integer                       :: m
  double precision              :: dem

! Output variables

  double precision,intent(out)   :: W(nBas,nBas,nBas,nBas)

!------------------------------------------------
! Compute static screening (physicist's notation)
!------------------------------------------------

  do p=1,nBas
    do q=1,nBas
      do r=1,nBas
        do s=1,nBas

          chi = 0d0
          do m=1,nS
            dem = Om(m)**2 + eta**2
            chi = chi + rho(p,q,m)*rho(r,s,m)*Om(m)/dem
          enddo

          W(p,s,q,r) = - lambda*ERI(p,s,q,r) + 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine 
