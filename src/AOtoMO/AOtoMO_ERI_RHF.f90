
! ---

subroutine AOtoMO_ERI_RHF(nBas_AOs, nBas_MOs, c, ERI_AO, ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  double precision,intent(in)   :: ERI_AO(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: c(nBas_AOs,nBas_MOs)

! Local variables

  double precision,allocatable  :: a1(:,:,:,:)
  double precision,allocatable  :: a2(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas_MOs,nBas_MOs,nBas_MOs,nBas_MOs)

! Memory allocation

  allocate(a2(nBas_AOs,nBas_AOs,nBas_AOs,nBas_MOs))
  allocate(a1(nBas_AOs,nBas_AOs,nBas_MOs,nBas_MOs))

! Four-index transform via semi-direct O(N^5) algorithm

  call dgemm( 'T', 'N', nBas_AOs*nBas_AOs*nBas_AOs, nBas_MOs, nBas_AOs, 1.d0 &
            , ERI_AO(1,1,1,1), nBas_AOs, c(1,1), nBas_AOs                    &
            , 0.d0, a2(1,1,1,1), nBas_AOs*nBas_AOs*nBas_AOs)

  call dgemm( 'T', 'N', nBas_AOs*nBas_AOs*nBas_MOs, nBas_MOs, nBas_AOs, 1.d0 &
            , a2(1,1,1,1), nBas_AOs, c(1,1), nBas_AOs                        &
            , 0.d0, a1(1,1,1,1), nBas_AOs*nBas_AOs*nBas_MOs)

  deallocate(a2)
  allocate(a2(nBas_AOs,nBas_MOs,nBas_MOs,nBas_MOs))

  call dgemm( 'T', 'N', nBas_AOs*nBas_MOs*nBas_MOs, nBas_MOs, nBas_AOs, 1.d0 &
            , a1(1,1,1,1), nBas_AOs, c(1,1), nBas_AOs                        &
            , 0.d0, a2(1,1,1,1), nBas_AOs*nBas_MOs*nBas_MOs)

  deallocate(a1)

  call dgemm( 'T', 'N', nBas_MOs*nBas_MOs*nBas_MOs, nBas_MOs, nBas_AOs, 1.d0 &
            , a2(1,1,1,1), nBas_AOs, c(1,1), nBas_AOs                        &
            , 0.d0, ERI_MO(1,1,1,1), nBas_MOs*nBas_MOs*nBas_MOs)

  deallocate(a2)

end subroutine 




