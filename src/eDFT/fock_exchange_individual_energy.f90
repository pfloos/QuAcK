subroutine fock_exchange_individual_energy(nBas,Pw,P,ERI,Ex)

! Compute the Fock exchange potential

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: Fx(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: Ex

! Compute HF exchange matrix

  allocate(Fx(nBas,nBas))

  call fock_exchange_potential(nBas,Pw(:,:),ERI(:,:,:,:),Fx(:,:))
  Ex =       trace_matrix(nBas,matmul(P(:,:),Fx(:,:))) &
     - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:),Fx(:,:)))

end subroutine fock_exchange_individual_energy
