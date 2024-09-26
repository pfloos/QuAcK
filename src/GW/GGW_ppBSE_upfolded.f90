subroutine GGW_ppBSE_upfolded(nOrb,nC,nO,nV,nR,nS,ERI,rho,Om,eGW)

! Upfolded ppBSE@GW (TDA only!)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: eGW(nOrb)

! Local variables

  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: m,ij,kl,ijm
  integer                       :: ab,cd,abm
  integer,parameter             :: maxH = 100
  double precision              :: tmp,tmp1,tmp2,tmp3,tmp4

  integer                       :: n1h,n1p,n2h,n2p,n1h1p,n3h1p,n3p1h,n2h2p,nH
  double precision,external     :: Kronecker_delta

  integer,allocatable           :: order(:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Upfolded ppBSE@GW Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation by default!'
  write(*,*)

! Dimension of the supermatrix

  n1h = nO
  n1p = nV

  n2h = nO*(nO-1)/2
  n2p = nV*(nV-1)/2

  n1h1p = n1h*n1p

  n3h1p = n2h*n1h1p
  n3p1h = n2p*n1h1p

  nH = n2p + n3p1h 
! nH = n2h + n3h1p 

! Memory allocation

  allocate(order(nH),H(nH,nH),X(nH,nH),OmBSE(nH),Z(nH))

! Initialization

  H(:,:) = 0d0

!----------------------------------------!
!        Compute BSE supermatrix         !
!----------------------------------------!
!                                        !
!     |  D    -M1   -M1   -M2   -M2  |   ! 
!     |                              |   ! 
!     | +M1    E1    0     0     0   |   ! 
!     |                              |   ! 
! H = | +M1    0     E2    0     0   |   ! 
!     |                              |   ! 
!     | +M2    0     0     E3    0   |   ! 
!     |                              |   ! 
!     | +M2    0     0     0     E4  |   ! 
!                                        !
!----------------------------------------!

  !---------!
  ! Block C !
  !---------!

  ab = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
    ab = ab + 1
 
    cd = 0
    do c=nO+1,nOrb-nR
      do d=c+1,nOrb-nR
      cd = cd + 1
 
          H(ab,cd) = (eGW(a) + eGW(b))*Kronecker_delta(a,c)*Kronecker_delta(b,d) &
                     + (ERI(a,b,c,d) - ERI(a,b,d,c))
 
        end do
      end do
 
    end do
  end do

  !---------!
  ! Block D !
  !---------!

! ij = 0
! do i=nC+1,nO
!   do j=i+1,nO
!   ij = ij + 1
!
!   kl = 0
!   do k=nC+1,nO
!     do l=k+1,nO
!     kl = kl + 1
!
!         H(ij,kl) = - (eGW(i) + eGW(j))*Kronecker_delta(i,k)*Kronecker_delta(j,l) &
!                    + (ERI(i,j,k,l) - ERI(i,j,l,k))
!
!       end do
!     end do
!
!   end do
! end do

  !-----------------!
  ! Coupling Blocks !
  !-----------------!

  abm = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      do m=1,nS
        abm = abm + 1

        cd = 0
        do c=nO+1,nOrb-nR
          do d=c+1,nOrb-nR
            cd = cd + 1

            tmp1 = Kronecker_delta(b,d)*rho(a,c,m)
            tmp2 = Kronecker_delta(b,c)*rho(a,d,m)
            tmp3 = Kronecker_delta(a,d)*rho(b,c,m)
            tmp4 = Kronecker_delta(a,c)*rho(b,d,m)

            H(n2p+0*n3p1h+abm,cd             ) = tmp1 + tmp2
            H(cd             ,n2p+0*n3p1h+abm) = tmp3 + tmp4

!           H(n2h+1*n3h1p+ijm,kl             ) = +tmp4
!           H(kl             ,n2h+1*n3h1p+ijm) = +tmp2
!          
!           H(n2h+2*n3h1p+ijm,kl             ) = +tmp1
!           H(kl             ,n2h+2*n3h1p+ijm) = +tmp4

!           H(n2h+3*n3h1p+ijm,kl             ) = +tmp3
!           H(kl             ,n2h+3*n3h1p+ijm) = +tmp1

          end do
        end do

      end do
    end do
  end do

! ijm = 0
! do i=nC+1,nO
!   do j=i+1,nO
!     do m=1,nS
!       ijm = ijm + 1

!       kl = 0
!       do k=nC+1,nO
!         do l=k+1,nO
!           kl = kl + 1

!           tmp1 = Kronecker_delta(j,l)*rho(i,k,m)
!           tmp2 = Kronecker_delta(j,k)*rho(i,l,m)
!           tmp3 = Kronecker_delta(i,l)*rho(j,k,m)
!           tmp4 = Kronecker_delta(i,k)*rho(j,l,m)

!           H(n2h+0*n3h1p+ijm,kl             ) = tmp1 - tmp2
!           H(kl             ,n2h+0*n3h1p+ijm) = tmp3 - tmp4
!          
!           H(n2h+1*n3h1p+ijm,kl             ) = +tmp4
!           H(kl             ,n2h+1*n3h1p+ijm) = +tmp2
!          
!           H(n2h+2*n3h1p+ijm,kl             ) = +tmp1
!           H(kl             ,n2h+2*n3h1p+ijm) = +tmp4

!           H(n2h+3*n3h1p+ijm,kl             ) = +tmp3
!           H(kl             ,n2h+3*n3h1p+ijm) = +tmp1

!         end do
!       end do

!     end do
!   end do
! end do

  !------------!
  ! Block 3p1h !
  !------------!

  abm = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      do m=1,nS
        abm = abm + 1

        tmp = eGW(a) + eGW(b) + Om(m)

        H(n2p+0*n3p1h+abm,n2p+0*n3p1h+abm) = tmp

      end do
    end do
  end do

  !------------!
  ! Block 3h1p !
  !------------!

! ijm = 0
! do i=nC+1,nO
!   do j=i+1,nO
!     do m=1,nS
!       ijm = ijm + 1

!       tmp = - eGW(i) - eGW(j) + Om(m)

!       H(n2h+0*n3h1p+ijm,n2h+0*n3h1p+ijm) = tmp
!       H(n2h+1*n3h1p+ijm,n2h+1*n3h1p+ijm) = tmp
!       H(n2h+2*n3h1p+ijm,n2h+2*n3h1p+ijm) = tmp
!       H(n2h+3*n3h1p+ijm,n2h+3*n3h1p+ijm) = tmp

!     end do
!   end do
! end do

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

  call diagonalize_general_matrix(nH,H,OmBSE,X)

  do s=1,nH
    order(s) = s
  end do

  call quick_sort(OmBSE,order,nH)
  call set_order(X,order,nH,nH)

!-----------------!
! Compute weights !
!-----------------!

  Z(:) = 0d0
  do s=1,nH
    do ij=1,n2h
      Z(s) = Z(s) + X(ij,s)**2
    end do
  end do

!--------------!
! Dump results !
!--------------!

  write(*,*)'-------------------------------------------'
  write(*,*)'  Upfolded ppBSE excitation energies (eV)  '
  write(*,*)'-------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','OmBSE (eV)','|','Z','|'
  write(*,*)'-------------------------------------------'

  do s=1,min(nH,maxH)
    if(Z(s) > 1d-7) &
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',OmBSE(s)*HaToeV,'|',Z(s),'|'
  end do

  write(*,*)'-------------------------------------------'
  write(*,*)

end subroutine 
