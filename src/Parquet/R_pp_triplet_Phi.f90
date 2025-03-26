subroutine R_pp_triplet_Phi(nOrb,nC,nR,nOO,nVV,ee_trip_Om,ee_trip_rho,hh_trip_Om,hh_trip_rho,pp_trip_Phi)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nR,nOO,nVV
  double precision,intent(in)   :: ee_trip_Om(nVV)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_trip_Om(nOO)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n

! Output variables
  double precision,intent(out)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  pp_trip_Phi(:,:,:,:) = 0d0

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nVV
                 pp_trip_Phi(p,q,r,s) = pp_trip_Phi(p,q,r,s) &
                      + ee_trip_rho(p,q,n)*ee_trip_rho(r,s,n)/ee_trip_Om(n)            
              end do

              do n=1,nOO
                 pp_trip_Phi(p,q,r,s) = pp_trip_Phi(p,q,r,s) &
                      - hh_trip_rho(p,q,n)*hh_trip_rho(r,s,n)/hh_trip_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_pp_triplet_Phi
