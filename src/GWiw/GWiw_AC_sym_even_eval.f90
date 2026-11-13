complex(kind=8) function GWiw_AC_sym_even_eval(nblk,coef,w)

! Evaluate the symmetric-even continued fraction built by GWiw_AC_sym_even_coeff at the complex frequency w
! (I. Duchemin, beDeft, SymmetricEvenAnalyticContinuation_eval)

  implicit none

! Input variables

  integer,intent(in)            :: nblk
  double precision,intent(in)   :: coef(5,nblk)
  complex(kind=8),intent(in)    :: w

! Local variables

  integer                       :: k
  complex(kind=8)               :: z2,A0,A1,B0,B1,Sf,an,bn,A2,B2
  double precision              :: nrm

  z2 = w*w
  A0 = (1d0,0d0)
  A1 = (0d0,0d0)
  B0 = (0d0,0d0)
  B1 = (1d0,0d0)
  Sf = (1d0,0d0)

  do k=1,nblk
    an = coef(1,k)*Sf
    bn = coef(2,k) + coef(3,k)*(z2-coef(5,k))
    Sf = (z2-coef(5,k))**2 + coef(4,k)
    A2 = bn*A1 + an*A0
    B2 = bn*B1 + an*B0
    A0 = A1
    A1 = A2
    B0 = B1
    B1 = B2
    nrm = 1d0/max(abs(A0),abs(A1),abs(B0),abs(B1))
    A0 = A0*nrm
    A1 = A1*nrm
    B0 = B0*nrm
    B1 = B1*nrm
  end do

  GWiw_AC_sym_even_eval = A1/B1

end function
