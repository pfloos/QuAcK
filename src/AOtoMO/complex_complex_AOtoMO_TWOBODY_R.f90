subroutine complex_complex_AOtoMO_TWOBODY_R(nBas,nOrb,c,TWOBODY_AO,TWOBODY_MO)

! AO to MO transformation of two-electron integrals 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  complex*16,intent(in)         :: TWOBODY_AO(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: c(nBas,nOrb)

! Local variables

  complex*16,allocatable        :: a1(:,:,:,:)
  complex*16,allocatable        :: a2(:,:,:,:)
  complex*16,allocatable        :: complex_TWOBODY_AO(:,:,:,:)
  integer                       :: p, q, r, s, mu, nu, lam, sig
  complex*16                    :: val

! Output variables

  complex*16,intent(out)        :: TWOBODY_MO(nOrb,nOrb,nOrb,nOrb)

  if(nBas < 300) then
       
        ! Perform fast but memory consuming transformation via zgemm 
        ! Memory allocation
        
          allocate(a2(nBas,nBas,nBas,nOrb))
          allocate(complex_TWOBODY_AO(nBas,nBas,nBas,nBas))
          
          complex_TWOBODY_AO = (1d0,0d0)*TWOBODY_AO
        
        ! Four-index transform via semi-direct O(N^5) algorithm
        
          call zgemm( 'T', 'N', nBas*nBas*nBas, nOrb, nBas, 1.d0 &
                    , complex_TWOBODY_AO(1,1,1,1), nBas, c(1,1), nBas&
                    , 0.d0, a2(1,1,1,1), nBas*nBas*nBas)
          deallocate(complex_TWOBODY_AO)
          allocate(a1(nBas,nBas,nOrb,nOrb))
          call zgemm( 'T', 'N', nBas*nBas*nOrb, nOrb, nBas, 1.d0 &
                    , a2(1,1,1,1), nBas, c(1,1), nBas            &
                    , 0.d0, a1(1,1,1,1), nBas*nBas*nOrb)
        
          deallocate(a2)
          allocate(a2(nBas,nOrb,nOrb,nOrb))
        
          call zgemm( 'T', 'N', nBas*nOrb*nOrb, nOrb, nBas, 1.d0 &
                    , a1(1,1,1,1), nBas, c(1,1), nBas            &
                    , 0.d0, a2(1,1,1,1), nBas*nOrb*nOrb)
        
          deallocate(a1)
        
          call zgemm( 'T', 'N', nOrb*nOrb*nOrb, nOrb, nBas, 1.d0 &
                    , a2(1,1,1,1), nBas, c(1,1), nBas            &
                    , 0.d0, TWOBODY_MO(1,1,1,1), nOrb*nOrb*nOrb)
        
          deallocate(a2)
    else
          ! Perform naive for-sloops for memory efficient but slow transformation
          ! Init unused matrices 
          allocate(a1(0,0,0,0))
          allocate(a2(0,0,0,0))
          allocate(complex_TWOBODY_AO(0,0,0,0))
          ! Initialize output
          TWOBODY_MO = (0.d0, 0.d0)

          !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(p,q,r,s,mu,nu,lam,sig,val) COLLAPSE(4) SCHEDULE(dynamic)
          do p = 1, nOrb
            do q = 1, nOrb
              do r = 1, nOrb
                do s = 1, nOrb
                  val = (0.d0, 0.d0)
                  do mu = 1, nBas
                    do nu = 1, nBas
                      do lam = 1, nBas
                        do sig = 1, nBas
                          val = val + c(mu,p) * c(nu,q) * TWOBODY_AO(mu,nu,lam,sig) * c(lam,r) * c(sig,s)
                        end do
                      end do
                    end do
                  end do
                  TWOBODY_MO(p,q,r,s) = val
                end do
              end do
            end do
          end do
          !$OMP END PARALLEL DO 
          deallocate(a1,a2,complex_TWOBODY_AO)
     endif
end subroutine
