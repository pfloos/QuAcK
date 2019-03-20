function norm_coeff(alpha,a)

  implicit none

! Input variables

  double precision,intent(in)   :: alpha
  integer,intent(in)            :: a(3)

! local variable
  double precision              :: pi,dfa(3),dfac
  integer                       :: atot

! Output variable
  double precision norm_coeff

  pi = 4d0*atan(1d0)
  atot = a(1) + a(2) + a(3)

  dfa(1) = dfac(2*a(1))/(2d0**a(1)*dfac(a(1)))
  dfa(2) = dfac(2*a(2))/(2d0**a(2)*dfac(a(2)))
  dfa(3) = dfac(2*a(3))/(2d0**a(3)*dfac(a(3)))


  norm_coeff = (2d0*alpha/pi)**(3d0/2d0)*(4d0*alpha)**atot
  norm_coeff = norm_coeff/(dfa(1)*dfa(2)*dfa(3))
  norm_coeff = sqrt(norm_coeff)

end function norm_coeff
