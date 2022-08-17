subroutine static_screening_WB_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,ERI,Omega,rho,WB)

! Compute the VVOO block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: chi
  double precision              :: eps
  integer                       :: a,b,i,j,ab,ij,m

! Output variables

  double precision,intent(out)  :: WB(nVV,nOO)

! Initialization

  WB(:,:) = 0d0

!---------------!
! Singlet block !
!---------------!

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
            do m=1,nS
              eps = Omega(m)**2 + eta**2
              chi = chi + rho(a,j,m)*rho(b,i,m)*Omega(m)/eps
            enddo
           
            WB(ab,ij) = WB(ab,ij) + 4d0*lambda*chi/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(i,j)))

          end do
        end do
      end do
    end do

  end if

!---------------!
! Triplet block !
!---------------!

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
            do m=1,nS
              eps = Omega(m)**2 + eta**2
              chi = chi + rho(a,j,m)*rho(b,i,m)*Omega(m)/eps
            enddo

            WB(ab,ij) = WB(ab,ij) - 4d0*lambda*chi

          end do
        end do
      end do
    end do

  end if

end subroutine static_screening_WB_pp
