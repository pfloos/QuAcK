subroutine RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,nOOx,nVVx,lambda,Om1,rho1,Om2,rho2,TD)

! Compute the OOOO block of the static T-matrix 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  integer,intent(in)            :: nOOx
  integer,intent(in)            :: nVVx
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,k,l,ij,kl,ef,mn

! Output variables

  double precision,intent(out)  :: TD(nOOx,nOOx)

!===============!
! singlet block !
!===============!

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
      do j=i,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
 
            chi = 0d0
 
            do ef=1,nVV
              eps = + Om1(ef)
              chi = chi + rho1(i,j,ef)*rho1(k,l,ef)*eps/(eps**2 + eta**2)
            end do
 
            do mn=1,nOO
              eps = - Om2(mn)
              chi = chi + rho2(i,j,mn)*rho2(k,l,mn)*eps/(eps**2 + eta**2)
            end do
 
            TD(ij,kl) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

!===============!
! triplet block !
!===============!

  if(ispin == 2 .or. ispin == 4) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
 
            chi = 0d0
 
            do ef=1,nVV
              eps = + Om1(ef)
              chi = chi + rho1(i,j,ef)*rho1(k,l,ef)*eps/(eps**2 + eta**2)
            end do
 
            do mn=1,nOO
              eps = - Om2(mn)
              chi = chi + rho2(i,j,mn)*rho2(k,l,mn)*eps/(eps**2 + eta**2)
            end do
 
            TD(ij,kl) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

!==================!
! alpha-beta block !
!==================!

  if(ispin == 3) then

    ij = 0
    do i=nC+1,nO
      do j=nC+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=nC+1,nO
            kl = kl + 1
 
            chi = 0d0
 
            do ef=1,nVV
              eps = + Om1(ef)
              chi = chi + rho1(i,j,ef)*rho1(k,l,ef)*eps/(eps**2 + eta**2)
            end do
 
            do mn=1,nOO
              eps = - Om2(mn)
              chi = chi + rho2(i,j,mn)*rho2(k,l,mn)*eps/(eps**2 + eta**2)
            end do
 
            TD(ij,kl) = lambda*chi
 
          end do
        end do

      end do
    end do


  end if

end subroutine 
