subroutine G_linDM_GT(nOrb,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2,eta,linDM)
  
! Perform G0W0 calculation

implicit none
include 'parameters.h'


! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO,nVV
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om1(nVV),Om2(nOO)
  double precision,intent(in)   :: rho1(nOrb,nOrb,nVV)
  double precision,intent(in)   :: rho2(nOrb,nOrb,nOO)
  integer,intent(in)            :: eta

! Local variables

  integer                       :: i,j,k,a,b,c,n
  double precision              :: dem
  double precision              :: num
  
! Output variables

  double precision,intent(inout)  :: linDM(nOrb,nOrb)

  linDM(:,:) = 0d0

! OccOcc block of the density matrix
  
  do i=nC+1,nO
     do j=nC+1,nO
        do k=nC+1,nO
           do n=1,nVV

              num = - 1d0 *rho1(i,k,n)*rho1(j,k,n)
              dem = (Om1(n) - e(i) - e(k)) * (Om1(n) - e(j) - e(k))
              linDM(i,j) = linDM(i,j) + num/dem
              
           end do

        end do
     end do
  end do

  ! VirVir block of the density matrix
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        do c=nO+1,nOrb-nR
           do n=1,nOO

              num = rho2(a,c,n)*rho2(b,c,n)
              dem = (Om2(n) - e(a) - e(c)) * (Om2(n) - e(b) - e(c))
              linDM(a,b) = linDM(a,b) + num/dem
              
           end do
        end do
     end do
  end do

  ! OccVir block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        
        do j=nC+1,nO
           do n=1,nVV

              num = rho1(i,j,n)*rho1(a,j,n)
              dem = (e(i) - e(a)) * (e(i) + e(j) - Om1(n))
              linDM(i,a) = linDM(i,a) + num/dem
              
           end do
        end do

        do b=nO+1,nOrb-nR
           do n=1,nOO

              num = - rho2(i,b,n)*rho2(a,b,n)
              dem = (e(a) - e(i)) * (e(a) + e(b) - Om2(n))
              linDM(i,a) = linDM(i,a) + num/dem
              
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
  
end subroutine G_linDM_GT
