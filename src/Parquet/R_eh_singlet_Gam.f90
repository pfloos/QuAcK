subroutine R_eh_singlet_Gamma_A(nOrb,nC,nO,nR,nS,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi,eh_sing_Gam_A)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n

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
              
              ! do n=1,nS
              !    eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
              !         + 0.5d0 * eh_sing_rho(b,a,n)*eh_sing_rho(j,i,n)/eh_sing_Om(n) &
              !         + 0.5d0 * eh_sing_rho(a,b,n)*eh_sing_rho(i,j,n)/eh_sing_Om(n) &
              !         + 1.5d0 * eh_trip_rho(b,a,n)*eh_trip_rho(j,i,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(a,b,n)*eh_trip_rho(i,j,n)/eh_trip_Om(n)     
              ! end do

              ! do n=1,nVVs
              !    eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
              !         + ee_sing_rho(a,j,n)*ee_sing_rho(i,b,n)/ee_sing_Om(n)            
              ! end do

              ! do n=1,nOOs
              !    eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
              !         - hh_sing_rho(a,j,n)*hh_sing_rho(i,b,n)/hh_sing_Om(n)           
              ! end do

              ! do n=1,nVVt
              !    eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
              !         + 3d0 * ee_trip_rho(a,j,n)*ee_trip_rho(i,b,n)/ee_trip_Om(n)            
              ! end do

              ! do n=1,nOOt
              !    eh_sing_Gam_A(ia,jb) = eh_sing_Gam_A(ia,jb) &
              !         - 3d0 * hh_trip_rho(a,j,n)*hh_trip_rho(i,b,n)/hh_trip_Om(n)             
              ! end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_singlet_Gamma_A

subroutine R_eh_singlet_Gamma_B(nOrb,nC,nO,nR,nS,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi,eh_sing_Gam_B)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb
  integer                       :: n

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
              
              ! do n=1,nS
              !    eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
              !         + 0.5d0 * eh_sing_rho(j,a,n)*eh_sing_rho(b,i,n)/eh_sing_Om(n) &
              !         + 0.5d0 * eh_sing_rho(a,j,n)*eh_sing_rho(i,b,n)/eh_sing_Om(n) &
              !         + 1.5d0 * eh_trip_rho(j,a,n)*eh_trip_rho(b,i,n)/eh_trip_Om(n) &
              !         + 1.5d0 * eh_trip_rho(a,j,n)*eh_trip_rho(i,b,n)/eh_trip_Om(n)     
              ! end do

              ! do n=1,nVVs
              !    eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
              !         + ee_sing_rho(a,b,n)*ee_sing_rho(i,j,n)/ee_sing_Om(n)            
              ! end do

              ! do n=1,nOOs
              !    eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
              !         - hh_sing_rho(a,b,n)*hh_sing_rho(i,j,n)/hh_sing_Om(n)           
              ! end do

              ! do n=1,nVVt
              !    eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
              !         + 3d0 * ee_trip_rho(a,b,n)*ee_trip_rho(i,j,n)/ee_trip_Om(n)            
              ! end do

              ! do n=1,nOOt
              !    eh_sing_Gam_B(ia,jb) = eh_sing_Gam_B(ia,jb) &
              !         - 3d0 * hh_trip_rho(a,b,n)*hh_trip_rho(i,j,n)/hh_trip_Om(n)             
              ! end do
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_singlet_Gamma_B
