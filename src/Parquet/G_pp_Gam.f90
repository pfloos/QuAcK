subroutine G_pp_Gamma_D(nOrb,nC,nO,nOO,eh_Phi,pp_Gam_D)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nOO
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: ij,kl

! Output variables
  double precision, intent(out) :: pp_Gam_D(nOO,nOO)

! Initialization
  pp_Gam_D(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, i, j, ij, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_trip_Gam_B, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)
  ij = 0
  do i=nC+1,nO
     do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
           do l=k+1,nO
              kl = kl +1

              pp_Gam_D(ij,kl) = eh_Phi(i,j,k,l) - eh_Phi(i,j,l,k) 

           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine G_pp_Gamma_C(nOrb,nO,nR,nVV,eh_Phi,pp_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nVV
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: a0, aa
  integer                       :: ab,cd

! Output variables
  double precision, intent(out) :: pp_Gam_C(nVV,nVV)

! Initialization
  pp_Gam_C(:,:) = 0d0
  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)                   &
  !$OMP          PRIVATE(a, b, aa, ab, c, d, cd) &
  !$OMP          SHARED(nO, nOrb, nR, a0, eh_Phi, pp_Gam_C)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b

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

subroutine G_pp_Gamma_B(nOrb,nC,nO,nR,nOO,nVV,eh_Phi,pp_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nOO,nVV
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: a0, aa
  integer                       :: ab,ij

! Output variables
  double precision, intent(out) :: pp_Gam_B(nVV,nOO)

! Initialization
  pp_Gam_B(:,:) = 0d0
  a0 = nOrb - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)                   &
  !$OMP          PRIVATE(a, b, aa, ab, i, j, ij) &
  !$OMP          SHARED(nO, nC, nOrb, nR, a0, eh_Phi, pp_Gam_B)
  !$OMP DO
  do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
        ab = aa + b
        
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
