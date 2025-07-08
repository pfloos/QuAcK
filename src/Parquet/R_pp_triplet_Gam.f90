subroutine R_pp_triplet_Gamma_D(nOrb,nC,nO,nOOt,eh_sing_Phi,eh_trip_Phi,pp_trip_Gam_D)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nOOt
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: ij,kl
  integer                       :: n

! Output variables
  double precision, intent(out) :: pp_trip_Gam_D(nOOt,nOOt)

! Initialization
  pp_trip_Gam_D(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(i, j, ij, k, l, kl, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_trip_Gam_D, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ij = 0
  do i=nC+1,nO
     do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
           do l=k+1,nO
              kl = kl +1
              
              pp_trip_Gam_D(ij,kl) = 0.5d0*eh_sing_Phi(i,j,k,l) + 0.5d0*eh_trip_Phi(i,j,k,l) &
                                   - 0.5d0*eh_sing_Phi(i,j,l,k) - 0.5d0*eh_trip_Phi(i,j,l,k)

           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine 

subroutine R_pp_triplet_Gamma_C(nOrb,nO,nR,nVVt,eh_sing_Phi,eh_trip_Phi,pp_trip_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVVt
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd,a0,aa
  integer                       :: n

! Output variables
  double precision, intent(out) :: pp_trip_Gam_C(nVVt,nVVt)

! Initialization
  pp_trip_Gam_C(:,:) = 0d0

  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, c, d, cd) &
  !$OMP SHARED(nO, nOrb, nR, a0, pp_trip_Gam_C, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b

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

subroutine R_pp_triplet_Gamma_B(nOrb,nC,nO,nR,nOOt,nVVt,eh_sing_Phi,eh_trip_Phi,pp_trip_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOOt,nVVt
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij,a0,aa
  integer                       :: n

! Output variables
  double precision, intent(out) :: pp_trip_Gam_B(nVVt,nOOt)

! Initialization
  pp_trip_Gam_B(:,:) = 0d0

  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)          &
  !$OMP PRIVATE(a, b, aa, ab, i, j, ij) &
  !$OMP SHARED(nO, nC, nOrb, nR, a0, pp_trip_Gam_B, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b

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
