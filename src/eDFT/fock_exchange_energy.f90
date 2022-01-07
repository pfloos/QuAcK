subroutine fock_exchange_energy(nBas,P,Fx,Ex)

! Compute the (exact) Fock exchange energy

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: Fx(nBas,nBas)

! Local variables

  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: Ex

! Compute HF exchange energy

  Ex = 0.5d0*trace_matrix(nBas,matmul(P,Fx))

end subroutine fock_exchange_energy
