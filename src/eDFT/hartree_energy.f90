subroutine hartree_energy(nBas,P,J,EH)

! Compute the unrestricted version of the Hartree energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas,nspin)
  double precision,intent(in)   :: J(nBas,nBas,nspin)

! Local variables

  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EH(nsp)

! Compute the components of the Hartree energy

  EH(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
  EH(2) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2))) &
        + 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,1)))
  EH(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

end subroutine hartree_energy
