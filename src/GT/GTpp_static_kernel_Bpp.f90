subroutine GTpp_static_kernel_Bpp(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,nOOx,nVVx,lambda,Om1,rho1,Om2,rho2,TB)

! Compute the VVOO block of the static T-matrix

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
  integer                       :: i,j,a,b,ij,ab,cd,kl

! Output variables

  double precision,intent(out)  :: TB(nVVx,nOOx)

!===============!
! singlet block !
!===============!

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i,nO
            ij = ij + 1
 
            chi = 0d0
 
            do cd=1,nVV
              eps = + Om1(cd)
              chi = chi + rho1(a,b,cd)*rho1(i,j,cd)*eps/(eps**2 + eta**2)
            end do
 
            do kl=1,nOO
              eps = - Om2(kl)
              chi = chi + rho2(a,b,kl)*rho2(i,j,kl)*eps/(eps**2 + eta**2)
            end do
 
            TB(ab,ij) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

!===============!
! triplet block !
!===============!

  if(ispin == 2 .or. ispin == 4) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a+1,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=i+1,nO
            ij = ij + 1
 
            chi = 0d0
 
            do cd=1,nVV
              eps = + Om1(cd)
              chi = chi + rho1(a,b,cd)*rho1(i,j,cd)*eps/(eps**2 + eta**2)
            end do
 
            do kl=1,nOO
              eps = - Om2(kl)
              chi = chi + rho2(a,b,kl)*rho2(i,j,kl)*eps/(eps**2 + eta**2)
            end do
 
            TB(ab,ij) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

!==================!
! alpha-beta block !
!==================!

  if(ispin == 3) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR
        ab = ab + 1

        ij = 0
        do i=nC+1,nO
          do j=nC+1,nO
            ij = ij + 1
 
            chi = 0d0
 
            do cd=1,nVV
              eps = + Om1(cd)
              chi = chi + rho1(a,b,cd)*rho1(i,j,cd)*eps/(eps**2 + eta**2)
            end do
 
            do kl=1,nOO
              eps = - Om2(kl)
              chi = chi + rho2(a,b,kl)*rho2(i,j,kl)*eps/(eps**2 + eta**2)
            end do
 
            TB(ab,ij) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

end subroutine 
