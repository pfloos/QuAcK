subroutine linear_response_B_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,B_pp)

! Compute the B matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision,external     :: Kronecker_delta
  integer                       :: a,b,i,j,ab,ij

! Output variables

  double precision,intent(out)  :: B_pp(nVV,nOO)

! Build B matrix for the singlet manifold
 
  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=a,nBas-nR
        ab = ab + 1
        ij = 0
        do i=nC+1,nO
         do j=i,nO
            ij = ij + 1
 
            B_pp(ab,ij) = (ERI(a,b,i,j) + ERI(a,b,j,i))/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
 
          end do
        end do
      end do
    end do

  end if

! Build the B matrix for the triplet manifold

  if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=a+1,nBas-nR
        ab = ab + 1
        ij = 0
        do i=nC+1,nO
         do j=i+1,nO
            ij = ij + 1
 
            B_pp(ab,ij) = ERI(a,b,i,j) - ERI(a,b,j,i)
 
          end do
        end do
      end do
    end do

  end if

! Build the B matrix for the spinorbital basis

  if(ispin == 3) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=a+1,nBas-nR
        ab = ab + 1
        ij = 0
        do i=nC+1,nO
         do j=i+1,nO
            ij = ij + 1
 
            B_pp(ab,ij) = ERI(a,b,i,j) - ERI(a,b,j,i)
 
          end do
        end do
      end do
    end do

  end if

end subroutine linear_response_B_pp
