subroutine GGTpp_ppBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,Om1,rho1,Om2,rho2,T,KB_sta)

! Compute the VVOO of the T-matrix static pp kernel

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
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nOrb,nOrb,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nOrb,nOrb,nOO)
  double precision,intent(in)   :: T(nOrb,nOrb,nOrb,nOrb)

! Local variables

  double precision              :: dem,num
  integer                       :: p,q,r,s,e,m
  integer                       :: i,j,a,b
  integer                       :: ij,ab,kl,cd

! Output variables

  double precision,intent(out)  :: KB_sta(nVV,nOO)
  
! Initialization
  KB_sta(:,:) = 0d0

! Computing the kernel
! This is the same code as for the GF(2) kernel with elements T instead of ERI
  ab = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      ab = ab + 1

      ij = 0
      do i=nC+1,nO
        do j=i+1,nO
          ij = ij + 1

          do m=nC+1,nO
            do e=nO+1,nOrb-nR
   
               dem = eGF(m) - eGF(e)
               num =       T(a,m,i,e) * T(e,b,m,j) + T(a,e,i,m) * T(m,b,e,j)
               num = num - T(b,m,i,e) * T(e,a,m,j) - T(b,e,i,m) * T(m,a,e,j)
               KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)
              
            end do
          end do

        end do
      end do

    end do
  end do
  
end subroutine 
