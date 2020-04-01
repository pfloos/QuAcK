subroutine fock_exchange_individual_energy(nBas,Pw,P,ERI,Ex)

! Compute the Fock exchange potential

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  double precision,allocatable  :: Fx(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: Ex

! Compute HF exchange matrix

  allocate(Fx(nBas,nBas))

  Fx(:,:) = 0d0

  do si=1,nBas
    do la=1,nBas
      do nu=1,nBas
        do mu=1,nBas
          Fx(mu,nu) = Fx(mu,nu) + P(la,si)*ERI(mu,la,si,nu)
        enddo
      enddo
    enddo
  enddo

  Ex = -0.25d0*trace_matrix(nBas,matmul(P,Fx))

end subroutine fock_exchange_individual_energy
