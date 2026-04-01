subroutine RGW_phBSE_static_kernel_fullmat(nBas,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,W)

! Compute the second-order static BSE kernel for the resonant block (only for singlets!)

  implicit none
  include 'parameters.h'

! Input variables

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

! Output variables

  double precision,intent(out)   :: W(nBas,nBas,nBas,nBas)

!------------------------------------------------
! Compute static screening (physicist's notation)
!------------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do r=nC+1,nBas-nR
        do s=nC+1,nBas-nR

          chi = 0d0
          do m=1,nS
            chi = chi + rho(p,q,m)*rho(r,s,m)/Om(m)
          end do

          W(p,s,q,r) = 4d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
