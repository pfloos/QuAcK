subroutine RGW_phBSE_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,rho,KB)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  double precision              :: chi,eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: KB(nS,nS)

! Initialization

  KB(:,:) = 0d0

! Build dynamic B matrix

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
 
          chi = 0d0
          do kc=1,nS

            eps =     - OmRPA(kc) - (eGW(b) - eGW(i))
            chi = chi + rho(i,b,kc)*rho(j,a,kc)*eps/(eps**2 + eta**2)

            eps =     - OmRPA(kc) - (eGW(a) - eGW(j))
            chi = chi + rho(i,b,kc)*rho(j,a,kc)*eps/(eps**2 + eta**2)

          end do

          KB(ia,jb) = KB(ia,jb) - 2d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
