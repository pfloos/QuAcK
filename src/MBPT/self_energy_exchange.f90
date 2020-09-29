subroutine self_energy_exchange(nBas,c,P,G,SigmaX)

! Compute exchange part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas),P(nBas,nBas),G(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: SigmaX(nBas,nBas)

! Compute exchange part of the self-energy in the AO basis

  call exchange_matrix_AO_basis(nBas,P,G,SigmaX)

! Compute exchange part of the self-energy in the MO basis

  SigmaX = matmul(transpose(c),matmul(SigmaX,c))

end subroutine self_energy_exchange
