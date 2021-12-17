subroutine ufG0W0(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Unfold G0W0 equations

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

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: klc,kcd,ija,iab

  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: cGW(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**********************************************'
  write(*,*)'|        Unfolded G0W0 calculation           |'
  write(*,*)'**********************************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation for dynamic screening by default!'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nS
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),cGW(nH,nH),eGW(nH),Z(nH))

! Initialization

  H(:,:) = 0d0

!---------------------------!
!  Compute GW supermatrix   !
!---------------------------!
!                           !
!     |   F   V2h1p V2p1h | ! 
!     |                   | ! 
! H = | V2h1p C2h1p     0 | ! 
!     |                   | ! 
!     | V2p1h   0   C2p1h | ! 
!                           !
!---------------------------!

  !-------------!
  ! Block C2h1p !
  !-------------!

  ija = 0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR
        ija = ija + 1

        klc = 0
        do k=nC+1,nO
          do l=nC+1,nO
            do c=nO+1,nBas-nR
              klc = klc + 1

              H(1+ija,1+klc) & 
                = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
                - 2d0*ERI(j,c,a,l))*Kronecker_delta(i,k)

            end do
          end do
        end do

      end do
    end do
  end do

  !-------------!
  ! Block C2p1h !
  !-------------!

  iab = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR
        iab = iab + 1

        kcd = 0
        do k=nC+1,nO
          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR
              kcd = kcd + 1

              H(1+n2h1p+iab,1+n2h1p+kcd) &
                = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)

            end do
          end do
        end do

      end do
    end do
  end do

  do p=nC+1,nBas

  !---------!
  ! Block F !
  !---------!

    H(1,1) = eHF(p)

  !-------------!
  ! Block V2h1p !
  !-------------!

    klc = 0
    do k=nC+1,nO
      do l=nC+1,nO
        do c=nO+1,nBas-nR
          klc = klc + 1

          H(1       ,1+klc) = sqrt(2d0)*ERI(p,c,k,l)
          H(1+klc,1       ) = sqrt(2d0)*ERI(p,c,k,l)

        end do
      end do
    end do

  !-------------!
  ! Block V2p1h !
  !-------------!

    kcd = 0
    do k=nC+1,nO
      do c=nO+1,nBas-nR
        do d=nO+1,nBas-nR
          kcd = kcd + 1

          H(1              ,1+n2h1p+kcd) = sqrt(2d0)*ERI(p,k,d,c)
          H(1+n2h1p+kcd,1              ) = sqrt(2d0)*ERI(p,k,d,c)

        end do
      end do
    end do

  !-------------------------!
  ! Diagonalize supermatrix !
  !-------------------------!

   cGW(:,:) = H(:,:)
   call diagonalize_matrix(nH,cGW,eGW)

  !-----------------!
  ! Compute weights !
  !-----------------!

    Z(:) = 0d0
    do s=1,nH
      Z(s) = Z(s) + cGW(1,s)**2
    end do

  !--------------!
  ! Dump results !
  !--------------!

    write(*,*)'-------------------------------------------'
    write(*,'(A35,I3)')' G0W0 energies (eV) for orbital ',p
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP (eV)','|','Z','|'
    write(*,*)'-------------------------------------------'
 
    do s=1,nH
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
    enddo
 
    write(*,*)'-------------------------------------------'
    write(*,*)

  end do

end subroutine ufG0W0
