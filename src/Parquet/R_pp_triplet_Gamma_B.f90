subroutine R_pp_triplet_Gamma_B(nOrb,nC,nO,nR,nOOt,nVVt,eh_sing_Phi,eh_trip_Phi,pp_trip_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOOt,nVVt
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij,aa
  integer                       :: n

! Output variables
  double precision, intent(out) :: pp_trip_Gam_B(nVVt,nOOt)

! Initialization
  pp_trip_Gam_B(:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, i, j, ij) &
  !$OMP SHARED(nO, nC, nOrb, nR, pp_trip_Gam_B, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a=nO+1,nOrb-nR
     do b=a+1,nOrb-nR
        ab = (a-(nO+1))*(nOrb-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a

        ij = 0
        do i=nC+1,nO
           do j=i+1,nO
              ij = ij +1
              
              pp_trip_Gam_B(ab,ij) = 0.5d0*eh_sing_Phi(a,b,i,j) + 0.5d0*eh_trip_Phi(a,b,i,j) &
                                   - 0.5d0*eh_sing_Phi(a,b,j,i) - 0.5d0*eh_trip_Phi(a,b,j,i)
              
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
