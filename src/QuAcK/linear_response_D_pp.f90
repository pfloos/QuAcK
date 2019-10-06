subroutine linear_response_D_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,D_pp)

! Compute the D matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: eF
  double precision              :: delta_spin
  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,k,l,ij,kl

! Output variables

  double precision,intent(out)  :: D_pp(nOO,nOO)

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build D matrix

  eF = e(nO) + e(nO+1)
 
  ij = 0
  do i=nC+1,nO
   do j=i+1,nO
      ij = ij + 1
      kl = 0
      do k=nC+1,nO
       do l=k+1,nO
          kl = kl + 1

          D_pp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                        + (1d0 + delta_spin)*ERI(k,l,i,j) - (1d0 - delta_dRPA)*ERI(k,l,j,i)

        enddo
      enddo
    enddo
  enddo

end subroutine linear_response_D_pp
