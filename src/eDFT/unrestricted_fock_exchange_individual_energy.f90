subroutine unrestricted_fock_exchange_individual_energy(nBas,nEns,Pw,P,ERI,LZx,Ex)

! Compute the HF individual energy in the unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: Pw(nBas,nBas,nspin)
  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: Fx(:,:,:)
  double precision,external     :: trace_matrix

  integer                       :: iEns
  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: LZx(nspin)
  double precision,intent(out)  :: Ex(nspin,nEns)

! Compute HF exchange matrix

  allocate(Fx(nBas,nBas,nspin))

  do ispin=1,nspin

    call unrestricted_fock_exchange_potential(nBas,Pw(:,:,ispin),ERI,Fx(:,:,ispin))

    LZx(ispin) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,ispin),Fx(:,:,ispin)))

    do iEns=1,nEns
      Ex(ispin,iEns) = - 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin,iEns),Fx(:,:,ispin)))
    end do

  end do


end subroutine unrestricted_fock_exchange_individual_energy
