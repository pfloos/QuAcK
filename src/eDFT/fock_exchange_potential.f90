subroutine fock_exchange_potential(nBas,P,ERI,Fx)

! Compute the Fock exchange potential

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Compute HF exchange matrix

  Fx(:,:) = 0d0
  do si=1,nBas
    do la=1,nBas
      do nu=1,nBas
        do mu=1,nBas
          Fx(mu,nu) = Fx(mu,nu) - P(la,si)*ERI(mu,la,si,nu)
        enddo
      enddo
    enddo
  enddo

end subroutine fock_exchange_potential
