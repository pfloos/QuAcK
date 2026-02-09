subroutine G_pp_Gamma_B(nOrb,nC,nO,nR,nOO,nVV,eh_Phi,pp_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOO,nVV
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: aa
  integer                       :: ab,ij

! Output variables
  double precision, intent(out) :: pp_Gam_B(nVV,nOO)

! Initialization
  pp_Gam_B(:,:) = 0d0
  !$OMP PARALLEL DEFAULT(NONE)                   &
  !$OMP          PRIVATE(a, b, aa, ab, i, j, ij) &
  !$OMP          SHARED(nO, nC, nOrb, nR, eh_Phi, pp_Gam_B)
  !$OMP DO
  do a=nO+1,nOrb-nR
     do b=a+1,nOrb-nR
        ab = (a-(nO+1))*(nOrb-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a
        
        ij = 0
        do i=nC+1,nO
           do j=i+1,nO
              ij = ij + 1

              pp_Gam_B(ab,ij) = eh_Phi(a,b,i,j) - eh_Phi(a,b,j,i)
              
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
