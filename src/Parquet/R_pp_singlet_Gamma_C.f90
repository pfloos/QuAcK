subroutine R_pp_singlet_Gamma_C(nOrb,nO,nR,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd,aa,a0
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_sing_Gam_C(nVVs,nVVs)

! Initialization
  pp_sing_Gam_C(:,:) = 0d0

  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, c, d, cd) &
  !$OMP SHARED(nO, nOrb, nR, a0, pp_sing_Gam_C, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b
        
        cd = 0
        do c=nO+1,nOrb - nR
           do d=c,nOrb - nR
              cd = cd +1

              pp_sing_Gam_C(ab,cd) = 0.5d0*eh_sing_Phi(a,b,c,d) - 1.5d0*eh_trip_Phi(a,b,c,d) &
                                   + 0.5d0*eh_sing_Phi(a,b,d,c) - 1.5d0*eh_trip_Phi(a,b,d,c)

              pp_sing_Gam_C(ab,cd) = pp_sing_Gam_C(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
              
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
