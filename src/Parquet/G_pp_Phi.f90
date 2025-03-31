subroutine G_pp_Phi(nOrb,nC,nR,nOO,nVV,ee_Om,ee_rho,hh_Om,hh_rho,pp_Phi)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
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

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nVV
                 pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                      + 2d0 * ee_rho(p,q,n)*ee_rho(r,s,n)/ee_Om(n)            
              end do

              do n=1,nOO                 
                 pp_Phi(p,q,r,s) = pp_Phi(p,q,r,s) &
                      - 2d0 * hh_rho(p,q,n)*hh_rho(r,s,n)/hh_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_pp_Phi
