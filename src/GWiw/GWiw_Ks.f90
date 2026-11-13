complex(kind=8) function GWiw_Ks(a,b,s)

! Entire SRG kernel K_s(a,b) = (a+b) [1-exp(-(a^2+b^2) s)] / (a^2+b^2), evaluated at complex arguments

  implicit none

  complex(kind=8),intent(in)    :: a,b
  double precision,intent(in)   :: s

  complex(kind=8)               :: x,f

  x = (a*a+b*b)*s

  if(abs(x) < 1d-6) then
    f = 1d0 - x/2d0 + x*x/6d0
  else
    f = (1d0-exp(-x))/x
  end if

  GWiw_Ks = (a+b)*s*f

end function
