subroutine RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,eGF,Taaaa,Tabab,Tbaab,KB_sta)

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
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: Taaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tabab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tbaab(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ij,ab,cd,kl

! Output variables

  double precision,intent(out)  :: KB_sta(nVV,nOO)

! Initialization

  KB_sta(:,:) = 0d0
  
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
              eps = 0d0
              chi = chi + 0d0
            end do
 
            do kl=1,nOO
              eps = 0d0
              chi = chi + 0d0
            end do
 
            KB_sta(ab,ij) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

!===============!
! triplet block !
!===============!

  if(ispin == 2) then

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
              eps = 0d0
              chi = chi + 0d0
            end do
 
            do kl=1,nOO
              eps = 0d0
              chi = chi + 0d0
            end do
 
            KB_sta(ab,ij) = lambda*chi
 
          end do
        end do

      end do
    end do

  end if

end subroutine 
