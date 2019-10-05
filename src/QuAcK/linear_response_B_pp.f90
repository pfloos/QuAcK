subroutine linear_response_B_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,B_pp)

! Compute the B matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_spin
  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: a,b,i,j,ab,ij

! Output variables

  double precision,intent(out)  :: B_pp(nVV,nOO)

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build A matrix
 
  ab = 0
  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR
      ab = ab + 1
      ij = 0
      do i=nC+1,nO
        do j=nC+1,nO
          ij = ij + 1

          B_pp(ab,ij) = (1d0 + delta_spin)*ERI(a,b,i,j) - (1d0 - delta_dRPA)*ERI(a,b,j,j)

        enddo
      enddo
    enddo
  enddo

end subroutine linear_response_B_pp
