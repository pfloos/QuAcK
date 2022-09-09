subroutine static_screening_WD_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,ERI,Omega,rho,WD)

! Compute the OOOO block of the static screening W for the pp-BSE

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
  integer                       :: i,j,k,l,ij,kl,m

! Output variables

  double precision,intent(out)  :: WD(nOO,nOO)

! Initialization

  WD(:,:) = 0d0

!---------------!
! Singlet block !
!---------------!

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
            do m=1,nS
              eps = Omega(m)**2 + eta**2
              chi = chi + rho(i,k,m)*rho(j,l,m)*Omega(m)/eps &
                        + rho(i,l,m)*rho(j,k,m)*Omega(m)/eps
            enddo
 
            WD(ij,kl) = + 4d0*lambda*chi/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))

          end do
        end do
      end do
    end do

  end if

!---------------!
! Triplet block !
!---------------!

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1

            chi = 0d0
            do m=1,nS
              eps = Omega(m)**2 + eta**2
              chi = chi + rho(i,k,m)*rho(j,l,m)*Omega(m)/eps &
                        - rho(i,l,m)*rho(j,k,m)*Omega(m)/eps
            enddo
 
            WD(ij,kl) = + 4d0*lambda*chi

          end do
        end do
      end do
    end do

  end if

!---------------!
! SpinOrb block !
!---------------!

  if(ispin == 4) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1

            chi = 0d0
            do m=1,nS
              eps = Omega(m)**2 + eta**2
              chi = chi + rho(i,k,m)*rho(j,l,m)*Omega(m)/eps &
                        - rho(i,l,m)*rho(j,k,m)*Omega(m)/eps
            enddo
 
            WD(ij,kl) = + 2d0*lambda*chi

          end do
        end do
      end do
    end do

  end if

end subroutine static_screening_WD_pp
