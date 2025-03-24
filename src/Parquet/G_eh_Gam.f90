subroutine G_eh_Gamma(nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
           eh_Om,eh_rho,ee_Om,ee_rho,hh_Om,hh_rho, &
           eh_Gam)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_Gam(nOrb,nOrb,nOrb,nOrb)
  
! Initialization
  eh_Gam(:,:,:,:) = 0d0

  do s = nC+1, nOrb-nR
     do r = nC+1, nOrb-nR
        do q = nC+1, nOrb-nR
           do p = nC+1, nOrb-nR
              
              ! do n=1,nS
              !    eh_sing_Gam(p,q,r,s) = eh_sing_Gam(p,q,r,s) &
              !         + eh_sing_rho(s,p,n)*eh_sing_rho(q,r,n)/eh_sing_Om(n) &
              !         + 3d0 * eh_trip_rho(s,p,n)*eh_trip_rho(q,r,n)/eh_trip_Om(n)     
              ! end do

              ! do n=1,nVVs
              !    eh_sing_Gam(p,q,r,s) = eh_sing_Gam(p,q,r,s) &
              !         + ee_sing_rho(p,q,n)*ee_sing_rho(r,s,n)/ee_sing_Om(n)            
              ! end do

              ! do n=1,nOOs
              !    eh_sing_Gam(p,q,r,s) = eh_sing_Gam(p,q,r,s) &
              !         - hh_sing_rho(p,q,n)*hh_sing_rho(r,s,n)/hh_sing_Om(n)           
              ! end do

              ! do n=1,nVVt
              !    eh_sing_Gam(p,q,r,s) = eh_sing_Gam(p,q,r,s) &
              !         + 3d0 * ee_trip_rho(p,q,n)*ee_trip_rho(r,s,n)/ee_trip_Om(n)            
              ! end do

              ! do n=1,nOOt
              !    eh_sing_Gam(p,q,r,s) = eh_sing_Gam(p,q,r,s) &
              !         - 3d0 * hh_trip_rho(p,q,n)*hh_trip_rho(r,s,n)/hh_trip_Om(n)             
              ! end do
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_eh_Gamma

subroutine G_eh_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
           eh_Om,eh_rho,ee_Om,ee_rho,hh_Om,hh_rho, &
           eh_Gam_A)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_Gam_A(nS,nS)
  
! Initialization
  eh_Gam_A(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_Gam_A(ia,jb) = eh_Gam_A(ia,jb) &
                      +  eh_rho(b,a,n)*eh_rho(j,i,n)/eh_Om(n) &
                      +  eh_rho(a,b,n)*eh_rho(i,j,n)/eh_Om(n)     
              end do

              do n=1,nVV
                 eh_Gam_A(ia,jb) = eh_Gam_A(ia,jb) &
                      + 2d0 * ee_rho(a,j,n)*ee_rho(i,b,n)/ee_Om(n)            
              end do

              do n=1,nOO
                 eh_Gam_A(ia,jb) = eh_Gam_A(ia,jb) &
                      - 2d0 * hh_rho(a,j,n)*hh_rho(i,b,n)/hh_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_eh_Gamma_A

subroutine G_eh_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
           eh_Om,eh_rho,ee_Om,ee_rho,hh_Om,hh_rho, &
           eh_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n
  double precision,external     :: Kronecker_delta

! Output variables
  double precision, intent(out) :: eh_Gam_B(nS,nS)
  
! Initialization
  eh_Gam_B(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        ia = ia + 1
        
        jb = 0
        do j=nC+1,nO
           do b=nO+1,norb-nR
              jb = jb + 1
              
              do n=1,nS
                 eh_Gam_B(ia,jb) = eh_Gam_B(ia,jb) &
                      +  eh_rho(j,a,n)*eh_rho(b,i,n)/eh_Om(n) &
                      +  eh_rho(a,j,n)*eh_rho(i,b,n)/eh_Om(n)     
              end do

              do n=1,nVV
                 eh_Gam_B(ia,jb) = eh_Gam_B(ia,jb) &
                      + 2d0 * ee_rho(a,b,n)*ee_rho(i,j,n)/ee_Om(n)            
              end do

              do n=1,nOO
                 eh_Gam_B(ia,jb) = eh_Gam_B(ia,jb) &
                      - 2d0 * hh_rho(a,b,n)*hh_rho(i,j,n)/hh_Om(n)           
              end do
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_eh_Gamma_B
