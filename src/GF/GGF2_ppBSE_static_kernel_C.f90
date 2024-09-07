subroutine GGF2_ppBSE_static_kernel_C(eta,nBas,nC,nO,nV,nR,nVV,lambda,ERI,eGF,KC_sta)

! Compute the resonant part of the static BSE@GF2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d,e
  integer                       :: a0,aa,ab,cd
  
! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization

  KC_sta(:,:) = 0d0
  
! Second-order correlation kernel for the block C of the singlet manifold

  ab = 0
  do a=nO+1,nBas-nR
    do b=a+1,nBas-nR
      ab = ab + 1
 
      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1
 
          do m=nC+1,nO
            do e=nO+1,nBas-nR
  
              dem = eGF(m) - eGF(e)
              num =       (ERI(a,m,c,e) - ERI(a,m,e,c)) * (ERI(e,b,m,d) - ERI(e,b,d,m))
              num = num + (ERI(a,e,c,m) - ERI(a,e,m,c)) * (ERI(m,b,e,d) - ERI(m,b,d,e))
              num = num - (ERI(b,m,c,e) - ERI(b,m,e,c)) * (ERI(e,a,m,d) - ERI(e,a,d,m))
              num = num - (ERI(b,e,c,m) - ERI(b,e,m,c)) * (ERI(m,a,e,d) - ERI(m,a,d,e))
                                                                            
              KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
        
            end do
          end do
 
        end do
      end do
 
    end do
  end do

end subroutine 
