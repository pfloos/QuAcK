subroutine R_pp_singlet_Gamma_D(nOrb,nC,nO,nOOs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_D)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nOOs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: ij,kl
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_sing_Gam_D(nOOs,nOOs)

! Initialization
  pp_sing_Gam_D(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(i, j, ij, k, l, kl, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_sing_Gam_D, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ij = 0
  do i=nC+1,nO
     do j=i,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
           do l=k,nO
              kl = kl +1
              
              pp_sing_Gam_D(ij,kl) = 0.5d0*eh_sing_Phi(i,j,k,l) - 1.5d0*eh_trip_Phi(i,j,k,l) &
                                   - 1.5d0*eh_sing_Phi(i,j,l,k) + 0.5d0*eh_trip_Phi(i,j,l,k)

              pp_sing_Gam_D(ij,kl) = pp_sing_Gam_D(ij,kl)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine 

subroutine R_pp_singlet_Gamma_C(nOrb,nO,nR,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_sing_Gam_C(nVVs,nVVs)

! Initialization
  pp_sing_Gam_C(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, c, d, cd, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_sing_Gam_C, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ab = 0
  do a=nO+1,nOrb - nR
     do b=a,nOrb - nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nOrb - nR
           do d=c,nOrb - nR
              cd = cd +1
              

              pp_sing_Gam_C(ab,cd) = 0.5d0*eh_sing_Phi(a,b,c,d) - 1.5d0*eh_trip_Phi(a,b,c,d) &
                                   - 1.5d0*eh_sing_Phi(a,b,d,c) + 0.5d0*eh_trip_Phi(a,b,d,c)

              pp_sing_Gam_C(ab,cd) = pp_sing_Gam_C(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine R_pp_singlet_Gamma_B(nOrb,nC,nO,nR,nOOs,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOOs,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_sing_Gam_B(nVVs,nOOs)

! Initialization
  pp_sing_Gam_B(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, i, j, ij, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_sing_Gam_B, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ab = 0
  do a=nO+1,nOrb - nR
     do b=a,nOrb - nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
           do j=i,nO
              ij = ij +1
              
              pp_sing_Gam_B(ab,ij) = 0.5d0*eh_sing_Phi(a,b,i,j) - 1.5d0*eh_trip_Phi(a,b,i,j) &
                                   - 1.5d0*eh_sing_Phi(a,b,j,i) + 0.5d0*eh_trip_Phi(a,b,j,i)

              pp_sing_Gam_B(ab,ij) = pp_sing_Gam_B(ab,ij)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine
