subroutine unrestricted_hartree_individual_energy(nBas,nEns,Pw,P,ERI,LZH,EH)

! Compute the hartree contribution to the individual energies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: Pw(nBas,nBas,nspin)
  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: J(:,:,:)
  double precision,external     :: trace_matrix

  integer                       :: iEns
  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: LZH(nsp)
  double precision,intent(out)  :: EH(nsp,nEns)

! Compute HF exchange matrix

  allocate(J(nBas,nBas,nspin))

  do ispin=1,nspin
    call unrestricted_hartree_potential(nBas,Pw(:,:,ispin),ERI,J(:,:,ispin))
  end do

  LZH(1) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,1)))
  LZH(2) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,2))) &
           - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,1)))
  LZH(3) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,2)))

  do iEns=1,nEns

    EH(1,iEns) = trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,1)))
    EH(2,iEns) = trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,2))) &
               + trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,1)))
    EH(3,iEns) = trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,2)))

  end do

end subroutine unrestricted_hartree_individual_energy
