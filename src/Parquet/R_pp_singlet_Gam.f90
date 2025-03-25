subroutine R_pp_singlet_Gamma_D(nOrb,nC,nO,nS,nOOs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_D)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nS,nOOs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: ij,kl
  integer                       :: n

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
              
              ! do n=1,nS
              !    pp_sing_Gam_D(ij,kl) = pp_sing_Gam_D(ij,kl) &
              !         - 0.5d0 * eh_sing_rho(k,i,n)*eh_sing_rho(j,l,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(i,k,n)*eh_sing_rho(l,j,n)/eh_sing_Om(n) & 
              !         + 1.5d0 * eh_trip_rho(k,i,n)*eh_trip_rho(j,l,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(i,k,n)*eh_trip_rho(l,j,n)/eh_trip_Om(n) &
              !         - 0.5d0 * eh_sing_rho(l,i,n)*eh_sing_rho(j,k,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(i,l,n)*eh_sing_rho(k,j,n)/eh_sing_Om(n) &
              !         + 1.5d0 * eh_trip_rho(l,i,n)*eh_trip_rho(j,k,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(i,l,n)*eh_trip_rho(k,j,n)/eh_trip_Om(n)               
              ! end do

              ! pp_sing_Gam_D(ij,kl) = pp_sing_Gam_D(ij,kl)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine 

subroutine R_pp_singlet_Gamma_C(nOrb,nO,nR,nS,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nO,nR,nS,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd
  integer                       :: n

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
              
              ! do n=1,nS
              !    pp_sing_Gam_C(ab,cd) = pp_sing_Gam_C(ab,cd) &
              !         - 0.5d0 * eh_sing_rho(c,a,n)*eh_sing_rho(b,d,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(a,c,n)*eh_sing_rho(d,b,n)/eh_sing_Om(n) & 
              !         + 1.5d0 * eh_trip_rho(c,a,n)*eh_trip_rho(b,d,n)/eh_trip_Om(n) & 
              !         + 1.5d0 * eh_trip_rho(a,c,n)*eh_trip_rho(d,b,n)/eh_trip_Om(n) &
              !         - 0.5d0 * eh_sing_rho(d,a,n)*eh_sing_rho(b,c,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(a,d,n)*eh_sing_rho(c,b,n)/eh_sing_Om(n) &
              !         + 1.5d0 * eh_trip_rho(d,a,n)*eh_trip_rho(b,c,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(a,d,n)*eh_trip_rho(c,b,n)/eh_trip_Om(n)              
              ! end do

              ! pp_sing_Gam_C(ab,cd) = pp_sing_Gam_C(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine R_pp_singlet_Gamma_B(nOrb,nC,nO,nR,nS,nOOs,nVVs,eh_sing_Phi,eh_trip_Phi,pp_sing_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS,nOOs,nVVs
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij
  integer                       :: n

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
              
              ! do n=1,nS
              !    pp_sing_Gam_B(ab,ij) = pp_sing_Gam_B(ab,ij) &
              !         - 0.5d0 * eh_sing_rho(i,a,n)*eh_sing_rho(b,j,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(a,i,n)*eh_sing_rho(j,b,n)/eh_sing_Om(n) & 
              !         + 1.5d0 * eh_trip_rho(i,a,n)*eh_trip_rho(b,j,n)/eh_trip_Om(n) & 
              !         + 1.5d0 * eh_trip_rho(a,i,n)*eh_trip_rho(j,b,n)/eh_trip_Om(n) &
              !         - 0.5d0 * eh_sing_rho(j,a,n)*eh_sing_rho(b,i,n)/eh_sing_Om(n) &
              !         - 0.5d0 * eh_sing_rho(a,j,n)*eh_sing_rho(i,b,n)/eh_sing_Om(n) &
              !         + 1.5d0 * eh_trip_rho(j,a,n)*eh_trip_rho(b,i,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(a,j,n)*eh_trip_rho(i,b,n)/eh_trip_Om(n)                              
              ! end do

              ! pp_sing_Gam_B(ab,ij) = pp_sing_Gam_B(ab,ij)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine
