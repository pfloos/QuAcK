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
