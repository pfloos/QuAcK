subroutine R_pp_triplet_Gamma_C(nOrb,nO,nR,nVVt,eh_sing_Phi,eh_trip_Phi,pp_trip_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVVt
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd,aa
  integer                       :: n

! Output variables
  double precision, intent(out) :: pp_trip_Gam_C(nVVt,nVVt)

! Initialization
  pp_trip_Gam_C(:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, c, d, cd) &
  !$OMP SHARED(nO, nOrb, nR, pp_trip_Gam_C, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a=nO+1,nOrb-nR
     do b=a+1,nOrb-nR
        ab = (a-(nO+1))*(nOrb-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a


        cd = 0
        do c=nO+1,nOrb - nR
           do d=c+1,nOrb - nR
              cd = cd +1

              pp_trip_Gam_C(ab,cd) = 0.5d0*eh_sing_Phi(a,b,c,d) + 0.5d0*eh_trip_Phi(a,b,c,d) &
                                   - 0.5d0*eh_sing_Phi(a,b,d,c) - 0.5d0*eh_trip_Phi(a,b,d,c)
              
           end do
        end do
        
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
