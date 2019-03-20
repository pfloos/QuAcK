subroutine CalcOmYuk(maxm,ExpG,ExpY,fG,NormYSq,Om)

! Compute the 0^m for the screened Coulomb operator: (00|f12/r12|00)^m

  implicit none

! Input variables

  integer,intent(in)            :: maxm
  double precision,intent(in)   :: ExpG,ExpY,fG,NormYSq

! Local variables

  integer                       :: m,k
  double precision              :: pi,t,dbinom
  double precision,allocatable  :: Fm(:)

! Output variables

  double precision,intent(inout):: Om(0:maxm)

  allocate(Fm(0:maxm))

  pi = 4d0*atan(1d0)

! Campute generalized Boys functions

  t = (ExpY - fG)*NormYSq
  call CalcBoysF(maxm,t,Fm)

! Compute (00|00)^m

  do m=0,maxm
    Om(m) = 0d0
    do k=0,m
      Om(m) = Om(m) + dbinom(m,k)*(ExpY/ExpG)**k*Fm(k)
    end do
    Om(m) = (2d0/sqrt(pi))*sqrt(ExpY)*(fG/ExpG)*exp(-fG*NormYSq)*Om(m)
  end do

  deallocate(Fm)

end subroutine CalcOmYuk
