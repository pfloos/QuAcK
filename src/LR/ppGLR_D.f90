subroutine ppGLR_D(nBas,nC,nO,nV,nR,nOO,lambda,e,ERI,Dpp)

! Compute the D matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: eF
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,k,l,ij,kl

! Output variables

  double precision,intent(out)  :: Dpp(nOO,nOO)

! Define the chemical potential

! eF = e(nO) + e(nO+1)
  eF = 0d0
 
! Build the D matrix for the spin-orbital basis 

  ij = 0
  do i=nC+1,nO
   do j=i+1,nO
      ij = ij + 1
      kl = 0
      do k=nC+1,nO
       do l=k+1,nO
          kl = kl + 1

          Dpp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                       + lambda*(ERI(i,j,k,l) - ERI(i,j,l,k))

        end do
      end do
    end do
  end do

end subroutine 
