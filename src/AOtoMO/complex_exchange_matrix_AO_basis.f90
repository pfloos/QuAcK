subroutine complex_exchange_matrix_AO_basis(nBas,P,ERI,K)

! Compute exchange matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  complex*16,intent(in)         :: P(nBas,nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  complex*16,intent(out)        :: K(nBas,nBas)

  K(:,:) = cmplx(0.d0,0.d0,kind=8)

  do nu=1,nBas
    do si=1,nBas
      do la=1,nBas
        do mu=1,nBas
          K(mu,nu) = K(mu,nu) - P(la,si)*ERI(mu,la,si,nu)
        end do
      end do
    end do
  end do

end subroutine 

! ---



