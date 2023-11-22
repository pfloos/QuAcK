subroutine AOtoMO_ERI_RHF(nBas,c,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c(nBas,nBas)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)
  integer                       :: mu,nu,la,si
  integer                       :: i,j,k,l

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation

  allocate(scr(nBas,nBas,nBas,nBas))

! Four-index transform via semi-direct O(N^5) algorithm

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,ERI_AO,nBas,c(1,1),size(c,1),0d0,scr,nBas**3)
  
  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,scr,nBas,c(1,1),size(c,1),0d0,ERI_MO,nBas**3)

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,ERI_MO,nBas,c(1,1),size(c,1),0d0,scr,nBas**3)

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,scr,nBas,c(1,1),size(c,1),0d0,ERI_MO,nBas**3)

end subroutine 
