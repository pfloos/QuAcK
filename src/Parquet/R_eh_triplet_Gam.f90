subroutine R_eh_triplet_Gamma(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                          eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho, &
                          ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho, &
                          hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho, eh_trip_Gam)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_trip_Om(nVVs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nVVs)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_trip_Gam(nOrb,nOrb,nOrb,nOrb)

! Initialization
  eh_trip_Gam(:,:,:,:) = 0d0

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              do n=1,nS
                 eh_trip_Gam(p,q,r,s) = eh_trip_Gam(p,q,r,s) &
                      + eh_sing_rho(s,p,n)*eh_sing_rho(q,r,n)/eh_sing_Om(n) &
                      - eh_trip_rho(s,p,n)*eh_trip_rho(q,r,n)/eh_trip_Om(n)     
              end do
              
              do n=1,nVVs
                 eh_trip_Gam(p,q,r,s) = eh_trip_Gam(p,q,r,s) &
                      - 0d0*ee_sing_rho(p,q,n) * ee_sing_rho(r,s,n)/ee_sing_Om(n)            
              end do

              do n=1,nOOs
                 eh_trip_Gam(p,q,r,s) = eh_trip_Gam(p,q,r,s) &
                      + 0d0*hh_sing_rho(p,q,n) * hh_sing_rho(r,s,n)/hh_sing_Om(n)           
              end do

              do n=1,nVVt
                 eh_trip_Gam(p,q,r,s) = eh_trip_Gam(p,q,r,s) &
                      + 0d0*ee_trip_rho(p,q,n) * ee_trip_rho(r,s,n)/ee_trip_Om(n)            
              end do

              do n=1,nOOt
                 eh_trip_Gam(p,q,r,s) = eh_trip_Gam(p,q,r,s) &
                      - 0d0*hh_trip_rho(p,q,n) * hh_trip_rho(r,s,n)/hh_trip_Om(n)             
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_triplet_Gamma

subroutine R_eh_triplet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                          eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho, &
                          ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho, &
                          hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho, eh_trip_Gam_A)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_trip_Om(nVVs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nVVs)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_trip_Gam_A(nS,nS)

! Initialization
  eh_trip_Gam_A(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_trip_Gam_A(ia,jb) = eh_trip_Gam_A(ia,jb) &
                      + eh_sing_rho(b,a,n)*eh_sing_rho(j,i,n)/eh_sing_Om(n) &
                      - eh_trip_rho(b,a,n)*eh_trip_rho(j,i,n)/eh_trip_Om(n)     
              end do

              do n=1,nVVs
                 eh_trip_Gam_A(ia,jb) = eh_trip_Gam_A(ia,jb) &
                      - ee_sing_rho(a,j,n)*ee_sing_rho(i,b,n)/ee_sing_Om(n)            
              end do

              do n=1,nOOs
                 eh_trip_Gam_A(ia,jb) = eh_trip_Gam_A(ia,jb) &
                      + hh_sing_rho(a,j,n)*hh_sing_rho(i,b,n)/hh_sing_Om(n)           
              end do

              do n=1,nVVt
                 eh_trip_Gam_A(ia,jb) = eh_trip_Gam_A(ia,jb) &
                      + ee_trip_rho(a,j,n)*ee_trip_rho(i,b,n)/ee_trip_Om(n)            
              end do

              do n=1,nOOt
                 eh_trip_Gam_A(ia,jb) = eh_trip_Gam_A(ia,jb) &
                      - hh_trip_rho(a,j,n)*hh_trip_rho(i,b,n)/hh_trip_Om(n)             
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_triplet_Gamma_A

subroutine R_eh_triplet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                          eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho, &
                          ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho, &
                          hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho, eh_trip_Gam_B)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_trip_Om(nVVs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nVVs)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_trip_Gam_B(nS,nS)

! Initialization
  eh_trip_Gam_B(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_trip_Gam_B(ia,jb) = eh_trip_Gam_B(ia,jb) &
                      + eh_sing_rho(j,a,n)*eh_sing_rho(b,i,n)/eh_sing_Om(n) &
                      - eh_trip_rho(j,a,n)*eh_trip_rho(b,i,n)/eh_trip_Om(n)     
              end do

              do n=1,nVVs
                 eh_trip_Gam_B(ia,jb) = eh_trip_Gam_B(ia,jb) &
                      - ee_sing_rho(a,b,n)*ee_sing_rho(i,j,n)/ee_sing_Om(n)            
              end do

              do n=1,nOOs
                 eh_trip_Gam_B(ia,jb) = eh_trip_Gam_B(ia,jb) &
                      + hh_sing_rho(a,b,n)*hh_sing_rho(i,j,n)/hh_sing_Om(n)           
              end do

              do n=1,nVVt
                 eh_trip_Gam_B(ia,jb) = eh_trip_Gam_B(ia,jb) &
                      + ee_trip_rho(a,b,n)*ee_trip_rho(i,j,n)/ee_trip_Om(n)            
              end do

              do n=1,nOOt
                 eh_trip_Gam_B(ia,jb) = eh_trip_Gam_B(ia,jb) &
                      - hh_trip_rho(a,b,n)*hh_trip_rho(i,j,n)/hh_trip_Om(n)             
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_triplet_Gamma_B
