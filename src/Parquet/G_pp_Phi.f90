subroutine G_pp_Phi(eta,nOrb,nC,nR,nOO,nVV,ee_Om,ee_rho,hh_Om,hh_rho,pp_Phi)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nR,nOO,nVV
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision              :: tmp
  integer                       :: nGrid,g
  double precision              :: wmin,wmax,dw
  double precision,allocatable  :: w(:)

! Output variables
  double precision, intent(out) :: pp_Phi(nOrb,nOrb,nOrb,nOrb)
  
! Initialization
  pp_Phi(:,:,:,:) = 0d0
  
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
  !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, pp_Phi, ee_rho, ee_Om, hh_rho, hh_Om, nGrid, w)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              

              do g=1,nGrid
                 do n=1,nVV
                    tmp = w(g) - ee_Om(n)
                    pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                         + (ee_rho(p,q,n)*ee_rho(r,s,n)/tmp) * (1d0 - exp(- 2d0 * eta * tmp * tmp))
                 end do
              
                 do n=1,nOO
                    tmp = w(g) - hh_Om(n)
                    pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                         - (hh_rho(p,q,n)*hh_rho(r,s,n)/tmp) * (1d0 - exp(- 2d0 * eta * tmp * tmp))           
                 end do
              end do
              
              pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) / nGrid
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(w)
  
end subroutine 
