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

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR

              do n=1,nVV
                 pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
                      + ee_sing_rho(p,q,n)*ee_sing_rho(r,s,n)/ee_sing_Om(n)            
              end do

              do n=1,nOO
                 pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
                      - hh_sing_rho(p,q,n)*hh_sing_rho(r,s,n)/hh_sing_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_pp_singlet_Phi
