subroutine AOtoMO_ERI(bra,ket,nBas,c,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm
! bra and ket are the spin of (bra|ket) = (bra bra|ket ket) = <bra ket|bra ket>

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: bra
  integer,intent(in)            :: ket
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas),c(nBas,nBas,nspin)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,k,l

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation

  allocate(scr(nBas,nBas,nBas,nBas))

! Four-index transform via semi-direct O(N^5) algorithm

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,ERI_AO,nBas,c(1,1,bra),size(c,1),0d0,scr,nBas**3)
  
  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,scr,nBas,c(1,1,ket),size(c,1),0d0,ERI_MO,nBas**3)

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,ERI_MO,nBas,c(1,1,bra),size(c,1),0d0,scr,nBas**3)

  call dgemm ('T','N',nBas**3,nBas,nBas,1d0,scr,nBas,c(1,1,ket),size(c,1),0d0,ERI_MO,nBas**3)

end subroutine 
