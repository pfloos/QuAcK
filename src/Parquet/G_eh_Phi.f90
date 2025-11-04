subroutine G_eh_Phi(eta,nOrb,nC,nR,nS,eh_Om,eh_rho,eh_Phi)

! Compute irreducible vertex in the eh channel
  implicit none

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nR,nS
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision              :: tmp
  integer                       :: nGrid,g
  double precision              :: wmin,wmax,dw
  double precision,allocatable  :: w(:)

! Output variables
  double precision, intent(out) :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  
! Initialization
  eh_Phi(:,:,:,:) = 0d0
  
! Minimum and maximum frequency values
  nGrid = 21
  wmin = - 25d0
  wmax = + 25d0
  dw = (wmax - wmin)/dble(nGrid)
  
  allocate(w(nGrid))
  
  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, r, s, n, g, tmp) &
  !$OMP SHARED(eta, nC, nOrb, nR, nS, eh_Phi, eh_rho, eh_Om, nGrid, w)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do g=1,nGrid
                 do n=1,nS
                 
                    tmp = w(g) - eh_Om(n)
                    eh_Phi(p,q,r,s) = eh_Phi(p,q,r,s)                                                           &
                                    + ( (eh_rho(p,r,n) * eh_rho(s,q,n) + eh_rho(r,p,n) * eh_rho(q,s,n)) / tmp ) &
                                    * (1d0 - exp(- 2d0 * eta * tmp * tmp))
                 end do
              end do
              
              eh_Phi(p,q,r,s) = eh_Phi(p,q,r,s) / nGrid
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(w)
  
end subroutine 
