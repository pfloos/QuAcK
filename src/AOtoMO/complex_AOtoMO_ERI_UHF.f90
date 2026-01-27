subroutine complex_AOtoMO_ERI_UHF(bra,ket,nBas,c,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: bra,ket,nBas
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: c(nBas,nBas,nspin)

! Local variables

  complex*16,allocatable        :: a1(:,:,:,:)
  complex*16,allocatable        :: a2(:,:,:,:)
  complex*16,allocatable        :: complex_ERI_AO(:,:,:,:)
  integer                       :: p, q, r, s, mu, nu, lam, sig
  complex*16                    :: val

! Output variables

  complex*16,intent(out)        :: ERI_MO(nBas,nBas,nBas,nBas)

  allocate(a2(nBas,nBas,nBas,nBas))
  allocate(complex_ERI_AO(nBas,nBas,nBas,nBas))
  
  complex_ERI_AO = (1d0,0d0)*ERI_AO
  
  ! Four-index transform via semi-direct O(N^5) algorithm
  
  call zgemm( 'T', 'N', nBas**3, nBas, nBas, (1.d0,0.d0) &
            , complex_ERI_AO(1,1,1,1), nBas, c(1,1,bra), nBas&
            , (0.d0,0.d0), a2(1,1,1,1), nBas*nBas*nBas)
  deallocate(complex_ERI_AO)
  allocate(a1(nBas,nBas,nBas,nBas))
  call zgemm( 'T', 'N', nBas*nBas*nBas, nBas, nBas, (1.d0,0.d0) &
            , a2(1,1,1,1), nBas, c(1,1,ket), nBas            &
            , (0.d0,0.d0), a1(1,1,1,1), nBas*nBas*nBas)
  
  deallocate(a2)
  allocate(a2(nBas,nBas,nBas,nBas))
  
  call zgemm( 'T', 'N', nBas**3, nBas, nBas, (1.d0,0.d0) &
            , a1(1,1,1,1), nBas, c(1,1,bra), nBas            &
            , (0.d0,0.d0), a2(1,1,1,1), nBas*nBas*nBas)
  
  deallocate(a1)
  
  call zgemm( 'T', 'N', nBas**3, nBas, nBas, (1.d0,0.d0) &
            , a2(1,1,1,1), nBas, c(1,1,ket), nBas            &
            , (0.d0,0.d0), ERI_MO(1,1,1,1), nBas*nBas*nBas)
  
  deallocate(a2)
end subroutine
