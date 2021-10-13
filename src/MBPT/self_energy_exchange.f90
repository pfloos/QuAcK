subroutine self_energy_exchange(nBas,c,P,ERI,SigX)

! Compute exchange part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas),P(nBas,nBas),ERI(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: SigX(nBas,nBas)

! Compute exchange part of the self-energy in the AO basis

  call exchange_matrix_AO_basis(nBas,P,ERI,SigX)

! Compute exchange part of the self-energy in the MO basis

  SigX = matmul(transpose(c),matmul(SigX,c))

end subroutine self_energy_exchange
