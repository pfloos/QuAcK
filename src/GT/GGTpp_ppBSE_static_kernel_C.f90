subroutine GGTpp_ppBSE_static_kernel_C(eta,nOrb,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,Om1,rho1,Om2,rho2,T,KC_sta)

! Compute the VVVV block of the T-matrix static pp kernels

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
  integer                       :: a,b,c,d,k,l
  integer                       :: ab,kl,cd

! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization
  KC_sta(:,:) = 0d0
  
! Computing the kernel
! This is the same code as for the GF(2) kernel with elements T instead of ERI
  ab = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      ab = ab + 1
 
      cd = 0
      do c=nO+1,nOrb-nR
        do d=c+1,nOrb-nR
          cd = cd + 1
 
          do m=nC+1,nO
            do e=nO+1,nOrb-nR
  
               dem = eGF(m) - eGF(e)
               num =       T(a,m,c,e) * T(e,b,m,d) + T(a,e,c,m) * T(m,b,e,d)
               num = num - T(b,m,c,e) * T(e,a,m,d) - T(b,e,c,m) * T(m,a,e,d)
               KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
        
            end do
          end do
 
        end do
      end do
 
    end do
  end do

end subroutine 
