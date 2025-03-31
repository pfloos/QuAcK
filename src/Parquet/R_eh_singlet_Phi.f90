subroutine R_eh_singlet_Phi(nOrb,nC,nR,nS,eh_sing_Om,eh_sing_rho,eh_sing_Phi)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nR,nS
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n

! Output variables
  double precision,intent(out)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  eh_sing_Phi(:,:,:,:) = 0d0

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nS
                 eh_sing_Phi(p,q,r,s) = eh_sing_Phi(p,q,r,s)                &
                      - eh_sing_rho(r,p,n)*eh_sing_rho(q,s,n)/eh_sing_Om(n) &
                      - eh_sing_rho(p,r,n)*eh_sing_rho(s,q,n)/eh_sing_Om(n)     
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine 
