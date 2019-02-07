subroutine CalcOmErf(maxm,ExpY,fG,NormYSq,Om)

! Compute the 0^m for the long-range Coulomb operator: (00|erf(r)/r|00)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: ExpY,fG,NormYSq

! Local variables

  integer                       :: m
  double precision              :: pi,t
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om (0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Campute generalized Boys functions

  t = fG*NormYSq
  call CalcBoysF(maxm,t,Fm)

! Compute (00|00)^m

  do m=0,maxm
    Om(m) = (2d0/sqrt(pi))*sqrt(fG)*(fG/ExpY)**m*Fm(m)
  enddo

  deallocate(Fm)

end subroutine CalcOmErf
