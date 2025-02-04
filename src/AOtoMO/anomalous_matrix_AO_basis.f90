subroutine anomalous_matrix_AO_basis(nBas,sigma,Pa,ERI,L)

! Compute anomalous L matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: Pa(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: L(nBas,nBas)

  L(:,:) = 0d0

  do nu=1,nBas
    do si=1,nBas
      do la=1,nBas
        do mu=1,nBas
          L(mu,nu) = L(mu,nu) + sigma*Pa(la,si)*ERI(la,si,mu,nu)
        end do
      end do
    end do
  end do

end subroutine 

! ---

