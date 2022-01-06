subroutine exchange_matrix_MO_basis(nBas,c,P,G,K)

! Compute exchange matrix in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas),P(nBas,nBas)
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: K(nBas,nBas)

! Compute Hartree Hamiltonian in the AO basis

  call exchange_matrix_AO_basis(nBas,P,G,K)

! Transform Coulomb matrix in the MO basis 

  K = matmul(transpose(c),matmul(K,c))

end subroutine exchange_matrix_MO_basis
