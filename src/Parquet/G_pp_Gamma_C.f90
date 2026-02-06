subroutine G_pp_Gamma_C(nOrb,nO,nR,nVV,eh_Phi,pp_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVV
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: aa
  integer                       :: ab,cd

! Output variables
  double precision, intent(out) :: pp_Gam_C(nVV,nVV)

! Initialization
  pp_Gam_C(:,:) = 0d0
  !$OMP PARALLEL DEFAULT(NONE)                   &
  !$OMP          PRIVATE(a, b, aa, ab, c, d, cd) &
  !$OMP          SHARED(nO, nOrb, nR, eh_Phi, pp_Gam_C)
  !$OMP DO
  do a=nO+1,nOrb-nR
     do b=a+1,nOrb-nR
        ab = (a-(nO+1))*(nOrb-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a

        cd = 0
        do c=nO+1,nOrb - nR
           do d=c+1,nOrb - nR
              cd = cd +1
              
              pp_Gam_C(ab,cd) = eh_Phi(a,b,c,d) - eh_Phi(a,b,d,c)
              
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
