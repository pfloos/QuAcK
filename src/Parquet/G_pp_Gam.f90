subroutine G_pp_Gamma(nOrb,nC,nO,nV,nR,nS,eh_Om,eh_rho,pp_Gam)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_Gam(nOrb,nOrb,nOrb,nOrb)

! Initialization
  pp_Gam(:,:,:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, i, j, ij, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_trip_Gam_B, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nS
                 pp_Gam(p,q,r,s) = pp_Gam(p,q,r,s) &
                      - eh_rho(r,p,n)*eh_rho(q,s,n)/eh_Om(n) &
                      - eh_rho(p,r,n)*eh_rho(s,q,n)/eh_Om(n) &
                      + eh_rho(s,p,n)*eh_rho(q,r,n)/eh_Om(n) &
                      + eh_rho(p,s,n)*eh_rho(r,q,n)/eh_Om(n)            
              end do
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine G_pp_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOO,eh_Om,eh_rho,pp_Gam_D)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: ij,kl
  integer                       :: n
  double precision,external     :: Kronecker_delta

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
              
              do n=1,nS
                 pp_Gam_D(ij,kl) = pp_Gam_D(ij,kl) &
                      - eh_rho(k,i,n)*eh_rho(j,l,n)/eh_Om(n) &
                      - eh_rho(i,k,n)*eh_rho(l,j,n)/eh_Om(n) &
                      + eh_rho(l,i,n)*eh_rho(j,k,n)/eh_Om(n) &
                      + eh_rho(i,l,n)*eh_rho(k,j,n)/eh_Om(n)            
              end do
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine G_pp_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVV,eh_Om,eh_rho,pp_Gam_C)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nVV
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: a,b,c,d
  integer                       :: ab,cd
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_Gam_C(nVV,nVV)

! Initialization
  pp_Gam_C(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, i, j, ij, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_trip_Gam_B, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ab = 0
  do a=nO+1,nOrb - nR
     do b=a+1,nOrb - nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nOrb - nR
           do d=c+1,nOrb - nR
              cd = cd +1
              
              do n=1,nS
                 
                 pp_Gam_C(ab,cd) = pp_Gam_C(ab,cd) &
                      - eh_rho(c,a,n)*eh_rho(b,d,n)/eh_Om(n) &
                      - eh_rho(a,c,n)*eh_rho(d,b,n)/eh_Om(n) &
                      + eh_rho(d,a,n)*eh_rho(b,c,n)/eh_Om(n) &
                      + eh_rho(a,d,n)*eh_rho(c,b,n)/eh_Om(n)

              end do
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine

subroutine G_pp_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOO,nVV,eh_Om,eh_rho,pp_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

! Local variables
  integer                       :: a,b,i,j
  integer                       :: ab,ij
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: pp_Gam_B(nVV,nOO)

! Initialization
  pp_Gam_B(:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)         &
!  !$OMP PRIVATE(a, b, ab, i, j, ij, n) &
!  !$OMP SHARED(nC, nOrb, nO, nS, pp_trip_Gam_B, eh_sing_rho, eh_sing_Om, eh_trip_rho, eh_trip_Om)
!  !$OMP DO COLLAPSE(2)

  ab = 0
  do a=nO+1,nOrb - nR
     do b=a+1,nOrb - nR
        ab = ab + 1
        
        ij = 0
        do i=nC+1,nO
           do j=i+1,nO
              ij = ij + 1

              do n=1,nS
                 pp_Gam_B(ab,ij) = pp_Gam_B(ab,ij) &
                      - eh_rho(i,a,n)*eh_rho(b,j,n)/eh_Om(n) &
                      - eh_rho(a,i,n)*eh_rho(j,b,n)/eh_Om(n) &
                      + eh_rho(j,a,n)*eh_rho(b,i,n)/eh_Om(n) &
                      + eh_rho(a,j,n)*eh_rho(i,b,n)/eh_Om(n)            
              end do
              
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine
