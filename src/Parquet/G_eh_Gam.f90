subroutine G_eh_Gamma_A(nOrb,nC,nO,nR,nS,eh_Phi,pp_Phi,eh_Gam_A)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Phi(nOrb,nOrb,nOrb,nOrb)

! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb

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
              
              eh_Gam_A(ia,jb) = - eh_Phi(a,j,b,i) + pp_Phi(a,j,i,b)
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_eh_Gamma_A

subroutine G_eh_Gamma_B(nOrb,nC,nO,nR,nS,eh_Phi,pp_Phi,eh_Gam_B)

! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Phi(nOrb,nOrb,nOrb,nOrb)
  
! Local variables
  integer                       :: i,a,j,b
  integer                       :: ia,jb

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
              
              eh_Gam_B(ia,jb) = - eh_Phi(a,b,j,i) + pp_Phi(a,b,i,j)
              
           enddo
        enddo
     enddo
  enddo
  
end subroutine G_eh_Gamma_B
