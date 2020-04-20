subroutine linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A_lr)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_spin,delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: A_lr(nS,nS)

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build A matrix

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          A_lr(ia,jb) = (e(a) - e(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                      + (1d0 + delta_spin)*lambda*ERI(i,b,a,j) &
                      - (1d0 - delta_dRPA)*lambda*ERI(i,b,j,a)

        enddo
      enddo
    enddo
  enddo

end subroutine linear_response_A_matrix
