subroutine R_eh_singlet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                          eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho, &
                          ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho, &
                          hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho, eh_sing_Gam_A)


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
  double precision, intent(out) :: eh_sing_Gam_A(nS,nS)

! Initialization
  eh_sing_Gam_A(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
                      + eh_sing_rho(a,b,n)*eh_sing_rho(j,i,n)/eh_sing_Om(n) &
                      + 3d0 * eh_trip_rho(a,b,n)*eh_trip_rho(j,i,n)/eh_trip_Om(n)     
              end do

              do n=1,nVVs
                 eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
                      + ee_sing_rho(a,j,n)*ee_sing_rho(i,b,n)/ee_sing_Om(n)            
              end do

              do n=1,nOOs
                 eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
                      - hh_sing_rho(a,j,n)*hh_sing_rho(i,b,n)/hh_sing_Om(n)           
              end do

              do n=1,nVVt
                 eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
                      + 3d0 * ee_trip_rho(a,j,n)*ee_trip_rho(i,b,n)/ee_trip_Om(n)            
              end do

              do n=1,nOOt
                 eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
                      - 3d0 * hh_trip_rho(a,j,n)*hh_trip_rho(i,b,n)/hh_trip_Om(n)             
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_singlet_Gamma_A

subroutine R_eh_singlet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                          eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho, &
                          ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho, &
                          hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho, eh_sing_Gam_B)


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
  double precision, intent(out) :: eh_sing_Gam_B(nS,nS)

! Initialization
  eh_sing_Gam_B(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
                      + eh_sing_rho(a,j,n)*eh_sing_rho(b,i,n)/eh_sing_Om(n) &
                      + 3d0 * eh_trip_rho(a,j,n)*eh_trip_rho(b,i,n)/eh_trip_Om(n)     
              end do

              do n=1,nVVs
                 eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
                      + ee_sing_rho(a,b,n)*ee_sing_rho(i,j,n)/ee_sing_Om(n)            
              end do

              do n=1,nOOs
                 eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
                      - hh_sing_rho(a,b,n)*hh_sing_rho(i,j,n)/hh_sing_Om(n)           
              end do

              do n=1,nVVt
                 eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
                      + 3d0 * ee_trip_rho(a,b,n)*ee_trip_rho(i,j,n)/ee_trip_Om(n)            
              end do

              do n=1,nOOt
                 eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
                      - 3d0 * hh_trip_rho(a,b,n)*hh_trip_rho(i,j,n)/hh_trip_Om(n)             
              end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_singlet_Gamma_B
