subroutine CalcOmERI(maxm,ExpY,NormYSq,Om)

! Compute the 0^m for ERIs: (00|00)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: ExpY,NormYSq

! Local variables

  integer                       :: m
  double precision              :: pi,t
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om (0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Campute generalized Boys functions

  t = ExpY*NormYSq
  call CalcBoysF(maxm,t,Fm)

! Compute (00|00)^m

  do m=0,maxm
    Om(m) = (2d0/sqrt(pi))*sqrt(ExpY)*Fm(m)
  enddo

  deallocate(Fm)

end subroutine CalcOmERI
