subroutine R_pp_singlet_Phi(nOrb,nC,nR,nOO,nVV,ee_sing_Om,ee_sing_rho,hh_sing_Om,hh_sing_rho,pp_sing_Phi)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nR,nOO,nVV
  double precision,intent(in)   :: ee_sing_Om(nVV)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_sing_Om(nOO)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n

! Output variables
  double precision,intent(out)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  pp_sing_Phi(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, r, s, n) &
  !$OMP SHARED(nC, nOrb, nR, nVV, nOO, pp_sing_Phi, ee_sing_rho, ee_sing_Om, hh_sing_rho, hh_sing_Om)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR

              do n=1,nVV
                 pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
                      + 2d0 * ee_sing_rho(p,q,n)*ee_sing_rho(r,s,n)/ee_sing_Om(n)            
              end do

              do n=1,nOO
                 pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
                      - 2d0 * hh_sing_rho(p,q,n)*hh_sing_rho(r,s,n)/hh_sing_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine 
