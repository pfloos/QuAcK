subroutine complex_Hartree_matrix_AO_basis(nBas,P,ERI,H)

! Compute Hartree matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  complex*16,intent(in)         :: P(nBas,nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
 
! Output variables

  complex*16,intent(out)  :: H(nBas,nBas)

  H(:,:) = cmplx(0.d0,0.d0,kind=8)

  do si=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do mu=1,nBas
           H(mu,nu) = H(mu,nu) + P(la,si)*ERI(mu,la,nu,si)
        end do
      end do
    end do
  end do

end subroutine 

! ---
