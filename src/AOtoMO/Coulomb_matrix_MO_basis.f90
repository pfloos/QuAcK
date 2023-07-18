subroutine Coulomb_matrix_MO_basis(nBas,c,P,G,J)

! Compute Coulomb matrix in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas),P(nBas,nBas)
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: J(nBas,nBas)

! Compute Hartree Hamiltonian in the AO basis

  call Coulomb_matrix_AO_basis(nBas,P,G,J)

! Transform Coulomb matrix in the MO basis 

  J = matmul(transpose(c),matmul(J,c))

end subroutine 
