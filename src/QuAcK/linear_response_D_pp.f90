subroutine linear_response_D_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,D_pp)

! Compute the D matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: eF
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,k,l,ij,kl

! Output variables

  double precision,intent(out)  :: D_pp(nOO,nOO)

! Define the chemical potential

  eF = e(nO) + e(nO+1)
 
! Build the D matrix for the singlet manifold

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
     do j=nC+1,i
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
         do l=nC+1,k
            kl = kl + 1
 
            D_pp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                          + ERI(k,l,i,j) - ERI(k,l,j,i)
 
          end do
        end do
      end do
    end do

  end if

! Build the D matrix for the triplet manifold

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
     do j=nC+1,i-1
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
         do l=nC+1,k-1
            kl = kl + 1
 
            D_pp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                          + (ERI(k,l,i,j) + ERI(k,l,j,i))/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
 
          end do
        end do
      end do
    end do

  end if

end subroutine linear_response_D_pp
