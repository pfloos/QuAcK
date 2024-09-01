subroutine huckel_guess(nBas, nOrb, S, Hc, X, c)

!  Hickel guess 

  implicit none

! Input variables

  integer,intent(in)            :: nBas, nOrb
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)

! Local variables

  integer                       :: mu,nu
  double precision              :: a

  double precision,allocatable  :: F(:,:)

! Output variables

  double precision,intent(out)  :: c(nBas,nOrb)

! Memory allocation

  allocate(F(nBas,nBas))

! Extended Huckel parameter

  a = 1.75d0

! GWH approximation

  do mu = 1, nBas
    F(mu,mu) = Hc(mu,mu)
    do nu = mu+1, nBas

      F(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
      F(nu,mu) = F(mu,nu)

    end do
  end do
  
  call core_guess(nBas, nOrb, F, X, c)

  deallocate(F)

end subroutine
