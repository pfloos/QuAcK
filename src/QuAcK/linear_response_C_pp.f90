subroutine linear_response_C_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,C_pp)

! Compute the C matrix of the pp channel

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
  double precision              :: eF
  double precision,external     :: Kronecker_delta

  integer                       :: a,b,c,d,ab,cd

! Output variables

  double precision,intent(out)  :: C_pp(nVV,nVV)

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build C matrix
 
  eF = e(nO) + e(nO+1)

  ab = 0
  do a=nO+1,nBas-nR
   do b=a+1,nBas-nR
      ab = ab + 1
      cd = 0
      do c=nO+1,nBas-nR
       do d=c+1,nBas-nR
          cd = cd + 1

          C_pp(ab,cd) = + (e(a) + e(b) - eF)*Kronecker_delta(a,c)*Kronecker_delta(b,d) & 
                        + (1d0 + delta_spin)*ERI(a,b,c,d) - (1d0 - delta_dRPA)*ERI(a,b,d,c)

        enddo
      enddo
    enddo
  enddo

end subroutine linear_response_C_pp
