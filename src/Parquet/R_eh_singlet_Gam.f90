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
              
              eh_sing_Gam_A(ia,jb) = - 0.5d0*eh_sing_Phi(a,j,b,i) - 1.5d0*eh_trip_Phi(a,j,b,i) &
                                     + 0.5d0*pp_sing_Phi(a,j,i,b) + 1.5d0*pp_trip_Phi(a,j,i,b)
              
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
              
              eh_sing_Gam_B(ia,jb) = - 0.5d0*eh_sing_Phi(a,b,j,i) - 1.5d0*eh_trip_Phi(a,b,j,i) &
                                     + 0.5d0*pp_sing_Phi(a,b,i,j) + 1.5d0*pp_trip_Phi(a,b,i,j)
              
           enddo
        enddo
     enddo
  enddo

end subroutine R_eh_singlet_Gamma_B
