subroutine AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  double precision,allocatable  :: a1(:,:,:,:)
  double precision,allocatable  :: a2(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nOrb,nOrb,nOrb,nOrb)

! Memory allocation

  allocate(a2(nBas,nBas,nBas,nOrb))
  allocate(a1(nBas,nBas,nOrb,nOrb))

! Four-index transform via semi-direct O(N^5) algorithm

  call dgemm( 'T', 'N', nBas*nBas*nBas, nOrb, nBas, 1.d0 &
            , ERI_AO(1,1,1,1), nBas, c(1,1), nBas        &
            , 0.d0, a2(1,1,1,1), nBas*nBas*nBas)


  call dgemm( 'T', 'N', nBas*nBas*nOrb, nOrb, nBas, 1.d0 &
            , a2(1,1,1,1), nBas, c(1,1), nBas            &
            , 0.d0, a1(1,1,1,1), nBas*nBas*nOrb)

  deallocate(a2)
  allocate(a2(nBas,nOrb,nOrb,nOrb))

  call dgemm( 'T', 'N', nBas*nOrb*nOrb, nOrb, nBas, 1.d0 &
            , a1(1,1,1,1), nBas, c(1,1), nBas            &
            , 0.d0, a2(1,1,1,1), nBas*nOrb*nOrb)

  deallocate(a1)

  call dgemm( 'T', 'N', nOrb*nOrb*nOrb, nOrb, nBas, 1.d0 &
            , a2(1,1,1,1), nBas, c(1,1), nBas            &
            , 0.d0, ERI_MO(1,1,1,1), nOrb*nOrb*nOrb)

  deallocate(a2)
  
end subroutine 
