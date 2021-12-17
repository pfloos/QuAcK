subroutine ufBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF,eGW)

! Unfold BSE@GW equations

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGW(nBas)

! Local variables

  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb,kc,iajb,kcld
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
  write(*,*)'|        Unfolded BSE@GW calculation          |'
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

!---------------------------!
!  Compute GW supermatrix   !
!---------------------------!
!                           !
!     |  A   -Ve   -Vh  |   ! 
!     |                 |   ! 
! H = | +Vh C2h2p   0   |   ! 
!     |                 |   ! 
!     | +Ve   0   C2p2h |   ! 
!                           !
!---------------------------!

  !---------!
  ! Block A !
  !---------!

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
    ia = ia + 1

    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
      jb = jb + 1

        H(ia,jb) = (eGW(a) - eGW(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                 + 2d0*ERI(i,b,a,j) - ERI(i,b,j,a)

        end do
      end do

    end do
  end do

  !----------------!
  ! Blocks Vp & Ve !
  !----------------!

  iajb=0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do j=nC+1,nO
        do b=nO+1,nBas-nR
        iajb = iajb + 1

        kc = 0
        do k=nC+1,nO
          do c=nO+1,nBas-nR
            kc = kc + 1

            tmp = sqrt(2d0)*Kronecker_delta(k,j)*ERI(b,a,c,i)
            H(n1h1p+iajb,kc              ) = +tmp
            H(kc        ,n1h1p+n2h2p+iajb) = -tmp
           
            tmp = sqrt(2d0)*Kronecker_delta(b,c)*ERI(a,k,i,j)
            H(n1h1p+n2h2p+iajb,kc        ) = +tmp
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
    do a=nO+1,nBas-nR
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          iajb = iajb + 1

          kcld = 0
          do k=nC+1,nO
            do c=nO+1,nBas-nR
              do l=nC+1,nO
                do d=nO+1,nBas-nR
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

  do s=1,nH
    if(Z(s) > 1d-7) &
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',Om(s)*HaToeV,'|',Z(s),'|'
  enddo

  write(*,*)'-------------------------------------------'
  write(*,*)

end subroutine ufBSE
