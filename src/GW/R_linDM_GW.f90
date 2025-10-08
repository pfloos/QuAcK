subroutine R_linDM_GW(nOrb,nC,nO,nV,nR,nS,e,Om,rho,eta,linDM)
  
! Perform G0W0 calculation

implicit none
include 'parameters.h'


! Input variables

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
  double precision              :: dem1,dem2
  double precision              :: num
  
! Output variables

  double precision,intent(inout)  :: linDM(nOrb,nOrb)

  linDM(:,:) = 0d0

! OccOcc block of the density matrix
  
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nOrb-nR
           do n=1,nS

              num = - 4d0 *rho(i,a,n)*rho(j,a,n)
              dem1 = e(i) - e(a) - Om(n)
              dem2 = e(j) - e(a) - Om(n)
              linDM(i,j) = linDM(i,j) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              
           end do
        end do
     end do
  end do

  ! VirVir block of the density matrix
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        do i=nC+1,nO
           do n=1,nS

              num = 4d0 *rho(i,a,n)*rho(i,b,n)
              dem1 = e(i) - e(a) - Om(n)
              dem2 = e(i) - e(b) - Om(n)
              linDM(a,b) = linDM(a,b) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              
           end do
        end do
     end do
  end do

  ! OccVir block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        
        do j=nC+1,nO
           do n=1,nS

              num = - 4d0 *rho(a,j,n)*rho(i,j,n)
              dem1 = e(i) - e(a)
              dem2 = e(j) - e(a) - Om(n)
              linDM(i,a) = linDM(i,a) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              
           end do
        end do

        do b=nO+1,nOrb-nR
           do n=1,nS

              num = 4d0 *rho(i,b,n)*rho(a,b,n)
              dem1 = e(i) - e(a)
              dem2 = e(i) - e(b) - Om(n)
              linDM(i,a) = linDM(i,a) + num*(dem1*dem2 - eta**2)/(dem1**2 + eta**2)/(dem2**2 + eta**2)
              
           end do
        end do
        
     end do
  end do

  ! VirOcc block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        
        linDM(a,i) = linDM(i,a)
        
     end do
  end do
  
end subroutine R_linDM_GW
