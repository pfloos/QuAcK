subroutine R_pp_singlet_Gamma_B(nOrb,nC,nO,nR,nOOs,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOOs,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij,aa,a0
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_sing_Gam_B(nVVs,nOOs)

! Initialization
  pp_sing_Gam_B(:,:) = 0d0

  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, i, j, ij) &
  !$OMP SHARED(nO, nC, nOrb, nR, a0, pp_sing_Gam_B, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b

        ij = 0
        do i=nC+1,nO
           do j=i,nO
              ij = ij +1
              
              pp_sing_Gam_B(ab,ij) = 0.5d0*eh_sing_Phi(a,b,i,j) - 1.5d0*eh_trip_Phi(a,b,i,j) &
                                   + 0.5d0*eh_sing_Phi(a,b,j,i) - 1.5d0*eh_trip_Phi(a,b,j,i)

              pp_sing_Gam_B(ab,ij) = pp_sing_Gam_B(ab,ij)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
              
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine
