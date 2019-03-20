subroutine CalcOmNuc(maxm,ExpPQi,NormPQSq,Om)

! Compute (0|V|0)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: ExpPQi,NormPQSq

! Local variables

  integer                       :: m
  double precision              :: pi,dm,t
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om (0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Campute generalized Boys functions

  t = NormPQSq/ExpPQi
  call CalcBoysF(maxm,t,Fm)

! Compute (00|00)^m

  do m=0,maxm
    dm =dble(m)
    Om(m) = (2d0/sqrt(pi))*(1d0/ExpPQi)**(dm+0.5d0)*Fm(m)
  end do

  deallocate(Fm)

end subroutine CalcOmNuc
