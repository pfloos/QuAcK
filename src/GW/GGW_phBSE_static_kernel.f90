subroutine GGW_phBSE_static_kernel(nOrb,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,W)

! Compute the BSE static kernel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  double precision              :: chi
  integer                       :: p,q,r,s,m

! Output variables

  double precision,intent(out)  :: W(nOrb,nOrb,nOrb,nOrb)

! Compute static kernel

  do p=nC+1,nOrb-nR
    do q=nC+1,nOrb-nR
      do r=nC+1,nOrb-nR
        do s=nC+1,nOrb-nR

          chi = 0d0
          do m=1,nS
            chi = chi + rho(p,q,m)*rho(r,s,m)/Om(m)
          end do

          W(p,s,q,r) = 2d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
