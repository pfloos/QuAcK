subroutine density_matrix(nBas,ON,c,P)

! Compute density matrix based on the occupation numbers

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ON(nBas),c(nBas,nBas)

! Local variables

  integer                       :: mu,nu,i

! Output variables

  double precision,intent(out)  :: P(nBas,nBas)

  P(:,:) = 0d0

  do i=1,nBas
    do nu=1,nBas
      do mu=1,nBas
        P(mu,nu) = P(mu,nu) + 2d0*ON(i)*c(mu,i)*c(nu,i)
      enddo
    enddo
  enddo

end subroutine 
