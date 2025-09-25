subroutine R_linDM_GT(nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,eta,linDM)
  
! Perform G0W0 calculation

implicit none
include 'parameters.h'


! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om1s(nVVs),Om2s(nOOs)
  double precision,intent(in)   :: Om1t(nVVt),Om2t(nOOt)
  double precision,intent(in)   :: rho1s(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: rho2s(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: rho1t(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: rho2t(nOrb,nOrb,nOOt)
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
           do n=1,nVVs

              num = - 1d0/2d0 *rho1s(i,k,n)*rho1s(j,k,n)
              dem = (Om1s(n) - e(i) - e(k)) * (Om1s(n) - e(j) - e(k))
              linDM(i,j) = linDM(i,j) + num/dem
              
           end do

           do n=1,nVVt

              num = - 3d0/2d0 *rho1t(i,k,n)*rho1t(j,k,n)
              dem = (Om1t(n) - e(i) - e(k)) * (Om1t(n) - e(j) - e(k))
              linDM(i,j) = linDM(i,j) + num/dem
              
           end do
        end do
     end do
  end do

  ! VirVir block of the density matrix
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        do c=nO+1,nOrb-nR
           do n=1,nOOs

              num = 1d0/2d0 *rho2s(a,c,n)*rho2s(b,c,n)
              dem = (Om2s(n) - e(a) - e(c)) * (Om2s(n) - e(b) - e(c))
              linDM(a,b) = linDM(a,b) + num/dem
              
           end do

           
           do n=1,nOOt

              num = 3d0/2d0 *rho2t(a,c,n)*rho2t(b,c,n)
              dem = (Om2t(n) - e(a) - e(c)) * (Om2t(n) - e(b) - e(c))
              linDM(a,b) = linDM(a,b) + num/dem
              
           end do
        end do
     end do
  end do

  ! OccVir block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR
        
        do j=nC+1,nO
           do n=1,nVVs

              num = 1d0/2d0 *rho1s(i,j,n)*rho1s(a,j,n)
              dem = (e(i) - e(a)) * (e(i) + e(j) - Om1s(n))
              linDM(i,a) = linDM(i,a) + num/dem
              
           end do

           do n=1,nVVt

              num = 3d0/2d0 *rho1t(i,j,n)*rho1t(a,j,n)
              dem = (e(i) - e(a)) * (e(i) + e(j) - Om1t(n))
              linDM(i,a) = linDM(i,a) + num/dem
              
           end do
        end do

        do b=nO+1,nOrb-nR
           do n=1,nOOs

              num = -1d0/2d0 *rho2s(i,b,n)*rho2s(a,b,n)
              dem = (e(a) - e(i)) * (e(a) + e(b) - Om2s(n))
              linDM(i,a) = linDM(i,a) + num/dem
              
           end do

           
           do n=1,nOOt

              num = -3d0/2d0 *rho2t(i,b,n)*rho2t(a,b,n)
              dem = (e(a) - e(i)) * (e(a) + e(b) - Om2t(n))
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
  
end subroutine R_linDM_GT
