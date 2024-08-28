subroutine huckel_guess(nBas_AOs, nBas_MOs, S, Hc, X, c)

!  Hickel guess 

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs)

! Local variables

  integer                       :: mu,nu
  double precision              :: a

  double precision,allocatable  :: F(:,:)

! Output variables

  double precision,intent(out)  :: c(nBas_AOs,nBas_MOs)

! Memory allocation

  allocate(F(nBas_AOs,nBas_AOs))

! Extended Huckel parameter

  a = 1.75d0

! GWH approximation

  do mu = 1, nBas_AOs
    F(mu,mu) = Hc(mu,mu)
    do nu = mu+1, nBas_AOs

      F(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
      F(nu,mu) = F(mu,nu)

    end do
  end do
  
  call core_guess(nBas_AOs, nBas_MOs, F, X, c)

  deallocate(F)

end subroutine
