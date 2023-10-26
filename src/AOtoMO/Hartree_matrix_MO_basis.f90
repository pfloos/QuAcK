subroutine Hartree_matrix_MO_basis(nBas,c,P,G,H)

! Compute Hartree matrix in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas),P(nBas,nBas)
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: H(nBas,nBas)

! Compute Hartree matrix in the AO basis

  call Hartree_matrix_AO_basis(nBas,P,G,H)

! Transform Hartree matrix in the MO basis 

  H = matmul(transpose(c),matmul(H,c))

end subroutine 
