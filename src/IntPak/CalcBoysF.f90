!module c_functions
!  use iso_c_binding
!  interface
!    function gsl_sf_gamma_inc_P(a,t) bind(C, name="gsl_sf_gamma_inc_P")
!      use iso_c_binding, only: c_double
!      real(kind=c_double), value :: a,t
!      real(kind=c_double) :: gsl_sf_gamma_inc_P
!    end function gsl_sf_gamma_inc_P
!  end interface
!end module

subroutine CalcBoysF(maxm,t,Fm)
!  use c_functions
! Comute the generalized Boys function Fm(t) using Slatec library

  implicit none

! Input variables

  double precision,intent(in)   :: t
  integer,intent(in)            :: maxm

! Local variables

  integer                       :: m
  double precision              :: dm
  double precision              :: dgami


! Output variables

  double precision,intent(inout):: Fm(0:maxm)

  if(t == 0d0) then
    do m=0,maxm
      dm = dble(m)
      Fm(m) = 1d0/(2d0*dm+1d0)
     end do
   else
     do m=0,maxm
       dm = dble(m)
!       Fm(m) = gamma(dm+0.5d0)*gsl_sf_gamma_inc_P(dm+0.5d0,t)/(2d0*t**(dm+0.5d0))
       Fm(m) = dgami(dm+0.5d0,t)/(2d0*t**(dm+0.5d0))
     end do
   end if

end subroutine CalcBoysF
