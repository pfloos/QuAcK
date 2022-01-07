subroutine mo_fock_exchange_potential(nBas,c,Fx,Vx)

! Compute the exchange potential in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: Fx(nBas,nBas)

! Local variables

  integer                       :: mu,nu
  integer                       :: p

! Output variables

  double precision,intent(out)  :: Vx(nBas)

! Compute Vx

  Vx(:) = 0d0
  do p=1,nBas
    do mu=1,nBas
      do nu=1,nBas
        Vx(p) = Vx(p) + c(mu,p)*Fx(mu,nu)*c(nu,p)
      end do
    end do
  end do

end subroutine mo_fock_exchange_potential
