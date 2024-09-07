subroutine ppGLR_B(nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,Bpp)

! Compute the B matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision,external     :: Kronecker_delta
  integer                       :: a,b,i,j,ab,ij

! Output variables

  double precision,intent(out)  :: Bpp(nVV,nOO)

! Build the B matrix for the spin-orbital basis

  ab = 0
  do a=nO+1,nBas-nR
   do b=a+1,nBas-nR
      ab = ab + 1
      ij = 0
      do i=nC+1,nO
       do j=i+1,nO
          ij = ij + 1
 
          Bpp(ab,ij) = lambda*(ERI(a,b,i,j) - ERI(a,b,j,i))
 
        end do
      end do
    end do
  end do

end subroutine 
