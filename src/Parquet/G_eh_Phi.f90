subroutine G_eh_Phi(nOrb,nC,nR,nS,eh_Om,eh_rho,eh_Phi)

! Compute irreducible vertex in the eh channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nR,nS
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n

! Output variables
  double precision, intent(out) :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  
! Initialization
  eh_Phi(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, r, s, n) &
  !$OMP SHARED(nC, nOrb, nR, nS, eh_Phi, eh_rho, eh_Om)
  !$OMP DO COLLAPSE(2)
  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nS                 
                 eh_Phi(p,q,r,s) = eh_Phi(p,q,r,s)           &
                      - eh_rho(r,p,n)*eh_rho(q,s,n)/eh_Om(n) &
                      - eh_rho(p,r,n)*eh_rho(s,q,n)/eh_Om(n)     
              end do
              
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine 
