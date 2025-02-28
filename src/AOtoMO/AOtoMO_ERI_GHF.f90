subroutine AOtoMO_ERI_GHF(nBas,nBas2,c1,c2,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm
! bra and ket are the spin of (bra1 bra2|ket1 ket2)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c1(nBas,nBas2)
  double precision,intent(in)   :: c2(nBas,nBas2)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas2,nBas2,nBas2,nBas2)

! Memory allocation

  allocate(scr(nBas2,nBas2,nBas2,nBas2))

  call dgemm('T','N',nBas**3,nBas2,nBas,1d0,ERI_AO,nBas,c1(1,1),nBas,0d0,scr,nBas**3)

  call dgemm('T','N',nBas**2*nBas2,nBas2,nBas,1d0,scr,nBas,c2(1,1),nBas,0d0,ERI_MO,nBas**2*nBas2)

  call dgemm('T','N',nBas*nBas2**2,nBas2,nBas,1d0,ERI_MO,nBas,c1(1,1),nBas,0d0,scr,nBas*nBas2**2)

  call dgemm('T','N',nBas2**3,nBas2,nBas,1d0,scr,nBas,c2(1,1),nBas,0d0,ERI_MO,nBas2**3)

  deallocate(scr)

end subroutine 
