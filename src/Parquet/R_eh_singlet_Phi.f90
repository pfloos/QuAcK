subroutine R_eh_singlet_Phi(eta,nOrb,nC,nR,nS,eh_sing_Om,eh_sing_rho,omega,eh_sing_Phi)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nR,nS
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: omega

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision              :: tmp
  integer                       :: nGrid,g
  double precision              :: wmin,wmax,dw
  double precision,allocatable  :: w(:)
  double precision              :: dem

! Output variables
  double precision,intent(out)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  eh_sing_Phi(:,:,:,:) = 0d0
  
! Minimum and maximum frequency values
! nGrid = 100
! wmin = - 10d0
! wmax = + 10d0
! dw = (wmax - wmin)/dble(nGrid)
  
! allocate(w(nGrid))
  
! do g=1,nGrid
!   w(g) = wmin + dble(g)*dw
! end do

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, r, s, n, g, dem) &
  !$OMP SHARED(eta, nC, nOrb, nR, nS, eh_sing_Phi, eh_sing_rho, eh_sing_Om, omega, nGrid, w)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR

!             do g=1,nGrid
                 do n=1,nS
                    
                    dem = omega - eh_sing_Om(n)
!                   dem = w(g) - eh_sing_Om(n)
                    eh_sing_Phi(p,q,r,s) = eh_sing_Phi(p,q,r,s) + (eh_sing_rho(p,r,n)*eh_sing_rho(s,q,n) / dem) &
                         * (1d0 - exp(- 2d0 * eta * dem * dem))
                    
                    dem = omega + eh_sing_Om(n)
!                   dem = w(g) + eh_sing_Om(n)
                    eh_sing_Phi(p,q,r,s) = eh_sing_Phi(p,q,r,s) - (eh_sing_rho(r,p,n)*eh_sing_rho(q,s,n) / dem) &
                         * (1d0 - exp(- 2d0 * eta * dem * dem))
                    
                 end do
!             end do
              
!             eh_sing_Phi(p,q,r,s) = eh_sing_Phi(p,q,r,s) / nGrid
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

! deallocate(w)

end subroutine 
