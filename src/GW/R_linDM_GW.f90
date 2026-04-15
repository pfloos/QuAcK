subroutine R_linDM_GW(flow,nOrb,nC,nO,nV,nR,nS,e,Om,rho,eta,DM)
  
! Compute the linearized GW density matrix

implicit none
include 'parameters.h'


! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  integer,intent(in)            :: eta

! Local variables

  integer                       :: i,j,a,b,n
  double precision              :: num
  double precision              :: dem1,dem2
  double precision              :: reg1,reg2
  double precision              :: s
  
! Output variables

  double precision,intent(inout)  :: DM(nOrb,nOrb)

! Initialization

  DM(:,:) = 0d0
  s = flow

! Occupied-occupied block of the density matrix
  
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nOrb-nR
           do n=1,nS

              num = - 4d0*rho(i,a,n)*rho(j,a,n)
              dem1 = e(i) - e(a) - Om(n)
              dem2 = e(j) - e(a) - Om(n)
 
              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

!             DM(i,j) = DM(i,j) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              DM(i,j) = DM(i,j) + num*reg1*reg2
              
           end do
        end do
     end do
  end do

  ! Virtual-virtual block of the density matrix
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        do i=nC+1,nO
           do n=1,nS

              num = 4d0*rho(i,a,n)*rho(i,b,n)
              dem1 = e(i) - e(a) - Om(n)
              dem2 = e(i) - e(b) - Om(n)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

!             DM(a,b) = DM(a,b) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              DM(a,b) = DM(a,b) + num*reg1*reg2
              
           end do
        end do
     end do
  end do

  ! Occupied-virtual and virtual-occupied block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        
        do j=nC+1,nO
           do n=1,nS

              num = - 4d0*rho(a,j,n)*rho(i,j,n)
              dem1 = e(i) - e(a)
              dem2 = e(j) - e(a) - Om(n)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

!             DM(i,a) = DM(i,a) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
!             DM(a,i) = DM(a,i) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              DM(i,a) = DM(i,a) + num*reg1*reg2
              DM(a,i) = DM(a,i) + num*reg1*reg2
              
           end do
        end do

        do b=nO+1,nOrb-nR
           do n=1,nS

              num = 4d0*rho(i,b,n)*rho(a,b,n)
              dem1 = e(i) - e(a)
              dem2 = e(i) - e(b) - Om(n)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

!             DM(i,a) = DM(i,a) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
!             DM(a,i) = DM(a,i) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              DM(i,a) = DM(i,a) + num*reg1*reg2
              DM(a,i) = DM(a,i) + num*reg1*reg2
              
           end do
        end do
        
     end do
  end do

! call matout(nOrb,nOrb,DM)
  
end subroutine
