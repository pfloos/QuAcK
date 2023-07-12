subroutine ppLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,Bpp)

! Compute the B matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
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
 
            Bpp(ab,ij) = lambda*(ERI(a,b,i,j) + ERI(a,b,j,i))/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
 
          end do
        end do
      end do
    end do

  end if

! Build the B matrix for the triplet manifold, or alpha-alpha, or in the spin-orbital basis

  if(ispin == 2 .or. ispin == 4) then

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

  end if

! Build the alpha-beta block of the B matrix

  if(ispin == 3) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=nO+1,nBas-nR
        ab = ab + 1
        ij = 0
        do i=nC+1,nO
         do j=nC+1,nO
            ij = ij + 1
 
            Bpp(ab,ij) = lambda*ERI(a,b,i,j)
 
          end do
        end do
      end do
    end do

  end if

end subroutine 
