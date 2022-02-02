subroutine huckel_guess(nBas,S,Hc,X,c)

!  Hickel guess 

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)

! Local variables

  integer                       :: mu,nu
  double precision              :: a

  double precision,allocatable  :: F(:,:)

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

! Memory allocation

  allocate(F(nBas,nBas))

! Extended Huckel parameter

  a = 1.75d0

! GWH approximation

  do mu=1,nBas
    F(mu,mu) = Hc(mu,mu)
    do nu=mu+1,nBas

      F(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
      F(nu,mu) = F(mu,nu)

    enddo
  enddo
  
  call core_guess(nBas,F,X,c)

end subroutine
