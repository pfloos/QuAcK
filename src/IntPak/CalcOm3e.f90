subroutine CalcOm3e(maxm,delta0,delta1,Y1,Y0,Om)

! Compute the 0^m for ERIs: (00|00)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: delta0,delta1,Y0,Y1

! Local variables

  integer                       :: m
  double precision              :: pi,t,OG
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om (0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Calculate OG

  OG = (pi**4/delta0)**(3d0/2d0)*exp(-Y0)

! Campute generalized Boys functions

  t = delta1/(delta1-delta0)*(Y1-Y0)
  call CalcBoysF(maxm,t,Fm)

! Compute (000|000)^m

  do m=0,maxm
    Om(m) = (2d0/sqrt(pi))*OG*sqrt(delta0/(delta1-delta0))*(delta1/(delta1-delta0))**m
    Om(m) = Om(m)*Fm(m)
  enddo

  deallocate(Fm)

end subroutine CalcOm3e
