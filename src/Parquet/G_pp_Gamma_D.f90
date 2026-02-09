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
