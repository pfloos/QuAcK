subroutine GGTpp_phBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nOO,nVV,nS,lambda,ERI,eGF,Om1,rho1,Om2,rho2,KA_sta)

! 

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nOrb,nOrb,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nOrb,nOrb,nOO)


! Local variables

  double precision              :: dem,num,chi,eps
  integer                       :: p,q,r,s,e,m
  integer                       :: a,b,c,d,i,j,k,l
  integer                       :: ab,kl,cd,ia,jb

! Output variables

  double precision,intent(out)  :: KA_sta(nS,nS)

! Initialization
  KA_sta(:,:) = 0d0
  
! Computing the kernel
  ia = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nOrb-nR
          jb = jb + 1

          chi = 0d0
          
          do cd=1,nVV
  
            eps = Om1(cd)
            chi = chi - rho1(i,j,cd)*rho1(a,b,cd)*eps/(eps**2 + eta**2)
        
          end do

          do kl=1,nOO

            eps = Om2(kl)
            chi = chi + rho2(i,j,kl)*rho2(a,b,kl)*eps/(eps**2 + eta**2)
             
          end do

          KA_sta(ia,jb) = lambda*chi
 
        end do
      end do
 
    end do
  end do

end subroutine 
