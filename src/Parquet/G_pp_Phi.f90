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

! Output variables
  double precision, intent(out) :: pp_Phi(nOrb,nOrb,nOrb,nOrb)
  
! Initialization
  pp_Phi(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, r, s, n) &
  !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, pp_Phi, ee_rho, ee_Om, hh_rho, hh_Om)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nVV
                 pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                      - (ee_rho(p,q,n)*ee_rho(r,s,n)/ee_Om(n)) * (1d0 - exp(- 2d0 * eta * ee_Om(n) * ee_Om(n)))
              end do

              do n=1,nOO                 
                 pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                      + (hh_rho(p,q,n)*hh_rho(r,s,n)/hh_Om(n)) * (1d0 - exp(- 2d0 * eta * hh_Om(n) * hh_Om(n)))           
              end do
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine 
