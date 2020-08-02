subroutine unrestricted_density_matrix(nBas,nEns,nO,c,P,occnum)

! Calculate density matrices

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: occnum(2,2,3)


! Local variables

  integer                       :: ispin
  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: P(nBas,nBas,nspin,nEns)

!-------------------------------------------------------
!------------------------ GOK-UKS ----------------------
!-------------------------------------------------------
  iEns = 1
  do ispin=1,nspin
    P(:,:,ispin,iEns) =  matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
  end do

  iEns = 2

  if(nO(1) > 1) then 
    P(:,:,1,iEns) = matmul(c(:,1:nO(1)-1,1),transpose(c(:,1:nO(1)-1,1)))
  else
    P(:,:,1,iEns)=0.d0 
  end if
  P(:,:,1,iEns) = P(:,:,1,iEns) + occnum(1,1,iEns)* matmul(c(:,nO(1):nO(1),1),transpose(c(:,nO(1):nO(1),1))) &
                + occnum(1,2,iEns)* matmul(c(:,nO(1)+1:nO(1)+1,1),transpose(c(:,nO(1)+1:nO(1)+1,1)))
 if(nO(2) > 1) then
    P(:,:,2,iEns) = matmul(c(:,1:nO(2)-1,2),transpose(c(:,1:nO(2)-1,2)))
  else
    P(:,:,2,iEns)=0.d0
  end if

  P(:,:,2,iEns) = P(:,:,2,iEns) + occnum(2,1,iEns)* matmul(c(:,1:nO(2),2),transpose(c(:,1:nO(2),2))) &
                + occnum(2,2,iEns)*matmul(c(:,1:nO(2)+1,2),transpose(c(:,1:nO(2)+1,2)))

 
  iEns = 3  

  if(nO(1) > 1) then
    P(:,:,1,iEns) = matmul(c(:,1:nO(1)-1,1),transpose(c(:,1:nO(1)-1,1)))
  else
    P(:,:,1,iEns)=0.d0
  end if
  P(:,:,1,iEns) = P(:,:,1,iEns) + occnum(1,1,iEns)* matmul(c(:,nO(1):nO(1),1),transpose(c(:,nO(1):nO(1),1))) &
                + occnum(1,2,iEns)* matmul(c(:,nO(1)+1:nO(1)+1,1),transpose(c(:,nO(1)+1:nO(1)+1,1)))

  if(nO(2) > 1) then
    P(:,:,2,iEns) = matmul(c(:,1:nO(2)-1,2),transpose(c(:,1:nO(2)-1,2)))
  else
    P(:,:,2,iEns)=0.d0
  end if
  P(:,:,2,iEns) = P(:,:,2,iEns) + occnum(2,1,iEns)* matmul(c(:,nO(2):nO(2),2),transpose(c(:,nO(2):nO(2),2))) &
                + occnum(2,2,iEns)*matmul(c(:,nO(2)+1:nO(2)+1,2),transpose(c(:,nO(2)+1:nO(2)+1,2)))

!-------------------------------------------------------------------
!--------------- For eDFT_UKS (fundamental gap)---------------------
!-------------------------------------------------------------------

! N-electron ground state

!  iEns = 1
!  do ispin=1,nspin
!    P(:,:,ispin,iEns) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
!  end do

! (N-1)-electron state: remove spin-up electrons

!  iEns = 2
!  P(:,:,2,iEns) = matmul(c(:,1:nO(2),2),transpose(c(:,1:nO(2),2)))  
!  if (nO(1) > 1) then
!    P(:,:,1,iEns) = matmul(c(:,1:nO(1)-1,1),transpose(c(:,1:nO(1)-1,1)))
!  else
!    P(:,:,1,iEns) = 0.d0
!  end if
! (N+1)-electron state: remove spin-down electrons

!  iEns = 3 
!  P(:,:,2,iEns) = matmul(c(:,1:nO(2)+1,2),transpose(c(:,1:nO(2)+1,2)))  
!  P(:,:,1,iEns) = matmul(c(:,1:nO(1),1),transpose(c(:,1:nO(1),1)))

end subroutine unrestricted_density_matrix
