subroutine ppLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,lambda,e,ERI,Dpp)

! Compute the D matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nOrb),ERI(nOrb,nOrb,nOrb,nOrb) 
  
! Local variables

  double precision              :: eF
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,k,l,ij,kl

! Output variables

  double precision,intent(out)  :: Dpp(nOO,nOO)

! Define the chemical potential

  eF = e(nO) + e(nO+1)
! eF = 0d0
 
! Build the D matrix for the singlet manifold

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
     do j=i,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
         do l=k,nO
            kl = kl + 1
 
            Dpp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                         + lambda*(ERI(i,j,k,l) + ERI(i,j,l,k))/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
 
          end do
        end do
      end do
    end do

  end if

! Build the D matrix for the triplet or alpha-alpha manifold

  if(ispin == 2 .or. ispin == 4) then

    ij = 0
    do i=nC+1,nO
     do j=i+1,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
         do l=k+1,nO
            kl = kl + 1
 
            Dpp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                         + lambda*(ERI(i,j,k,l) - ERI(i,j,l,k))
 
          end do
        end do
      end do
    end do

  end if

! Build the alpha-beta block of the D matrix

  if(ispin == 3) then

    ij = 0
    do i=nC+1,nO
     do j=nC+1,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
         do l=nC+1,nO
            kl = kl + 1
 
            Dpp(ij,kl) = - (e(i) + e(j) - eF)*Kronecker_delta(i,k)*Kronecker_delta(j,l) & 
                         + lambda*ERI(i,j,k,l)
 
          end do
        end do
      end do
    end do

  end if

end subroutine 
