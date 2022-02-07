subroutine self_energy_exchange_diag(nBas,c,P,ERI,SigX)

! Compute the diagonal elements of the exchange part of the self-energy

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

  double precision,intent(out)  :: SigX(nBas)

! Compute exchange part of the self-energy in the AO basis

  allocate(Fx(nBas,nBas))
  call exchange_matrix_AO_basis(nBas,P,ERI,Fx)

! Compute exchange part of the self-energy in the MO basis

  SigX(:) = 0d0
  do q=1,nBas
    do mu=1,nBas
      do nu=1,nBas
        SigX(q) = SigX(q) + c(mu,q)*Fx(mu,nu)*c(nu,q) 
      end do
    end do
  end do

  deallocate(Fx)

end subroutine self_energy_exchange_diag
