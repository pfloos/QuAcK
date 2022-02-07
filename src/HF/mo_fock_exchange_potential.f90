subroutine mo_fock_exchange_potential(nBas,c,P,ERI,Vx)

! Compute the exchange potential in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu
  integer                       :: q
  double precision,allocatable  :: Fx(:,:)

! Output variables

  double precision,intent(out)  :: Vx(nBas)

! Compute Vx

  allocate(Fx(nBas,nBas))
  call exchange_matrix_AO_basis(nBas,P,ERI,Fx)

  Vx(:) = 0d0
  do q=1,nBas
    do mu=1,nBas
      do nu=1,nBas
        Vx(q) = Vx(q) + c(mu,q)*Fx(mu,nu)*c(nu,q)
      end do
    end do
  end do

  deallocate(Fx)

end subroutine mo_fock_exchange_potential
