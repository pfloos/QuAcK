subroutine RGW_ppBSE_upfolded(nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF,eGW)

! Upfolded ppBSE@GW (TDA triplets only!)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eGW(nOrb)

! Local variables

  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb,kc,iajb,kcld
  integer,parameter             :: maxH = 20
  double precision              :: tmp

  integer                       :: n1h1p,n2h2p,nH
  double precision,external     :: Kronecker_delta

  integer,allocatable           :: order(:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**********************************************'
  write(*,*)'|       Unfolded ppBSE@GW calculation        |'
  write(*,*)'**********************************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation by default!'
  write(*,*)

! Dimension of the supermatrix

  n1h1p = nO*nV
  n2h2p = nO*nO*nV*nV
     nH = n1h1p + n2h2p + n2h2p

! Memory allocation

  allocate(order(nH),H(nH,nH),X(nH,nH),Om(nH),Z(nH))

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
  ! Block D !
  !---------!

  ij = 0
  do i=nC+1,nO
    do j=i+1,nO
    ij = ij + 1

    kl = 0
    do k=nC+1,nO
      do l=k+1,nO
      kl = kl + 1

          H(ij,kl) = - (eGW(i) + eGW(j))*Kronecker_delta(i,k)*Kronecker_delta(j,l) &
                   + (ERI(i,j,k,l) - ERI(i,j,l,k))

        end do
      end do

    end do
  end do

  !----------------!
  ! Blocks M1      !
  !----------------!

  ijm = 0
  do i=nC+1,nO
    do j=i+1,nO
      do m=1,nS
        ijm = ijm + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1

            tmp = sqrt(2d0)*Kronecker_delta(k,j)*M(i,k,m)
            H(n2h+ijm,kl     ) = +tmp
            H(kl     ,n2h+ijm) = -tmp
           
            tmp = sqrt(2d0)*Kronecker_delta(k,j)*M(i,k,m)
            H(n2h+1*n3h1p+ijm,kl     ) = +tmp
            H(kl     ,n2h+n3h1p+ijm) = -tmp
           
            tmp = sqrt(2d0)*Kronecker_delta(b,c)*M(j,k,m)
            H(n2h+2*n3h1p+iajb,kc        ) = +tmp
            H(kc              ,n1h1p+iajb) = -tmp

            tmp = sqrt(2d0)*Kronecker_delta(b,c)*M(j,k,m)
            H(n2h+3*n2h2p+iajb,kc        ) = +tmp
            H(kc              ,n1h1p+iajb) = -tmp

            end do
          end do

        end do
      end do
    end do
  end do

  !-------------!
  ! Block 2h2p !
  !-------------!

  iajb = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do j=nC+1,nO
        do b=nO+1,nOrb-nR
          iajb = iajb + 1

          kcld = 0
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do l=nC+1,nO
                do d=nO+1,nOrb-nR
                kcld = kcld + 1

                tmp = ((eHF(a) + eGW(b) - eHF(i) - eGW(j))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                    + 2d0*ERI(a,k,i,c))*Kronecker_delta(j,l)*Kronecker_delta(b,d)
                H(n1h1p      +iajb,n1h1p      +kcld) = tmp
                H(n1h1p+n2h2p+iajb,n1h1p+n2h2p+kcld) = tmp

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

! call matout(nH,nH,H)

  call diagonalize_general_matrix(nH,H,Om,X)

  do s=1,nH
    order(s) = s
  end do

  call quick_sort(Om,order,nH)
  call set_order(X,order,nH,nH)

!-----------------!
! Compute weights !
!-----------------!

  Z(:) = 0d0
  do s=1,nH
    do ia=1,n1h1p
      Z(s) = Z(s) + X(ia,s)**2
    end do
  end do

!--------------!
! Dump results !
!--------------!

  write(*,*)'-------------------------------------------'
  write(*,*)'   unfolded BSE excitation energies (eV)   '
  write(*,*)'-------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','Omega (eV)','|','Z','|'
  write(*,*)'-------------------------------------------'

  do s=1,min(nH,maxH)
    if(Z(s) > 1d-7) &
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',Om(s)*HaToeV,'|',Z(s),'|'
  end do

  write(*,*)'-------------------------------------------'
  write(*,*)

end subroutine 
