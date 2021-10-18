subroutine ufGW(eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Unfold GW equations

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
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
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: klc,kcd,ija,iab

  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGW(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**********************************************'
  write(*,*)'|          Unfolded GW calculation           |'
  write(*,*)'**********************************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation for dynamic screening by default!'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nBas*nO*nS
  n2p1h = nBas*nV*nS
  nH = nBas + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH))

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

  !---------!
  ! Block F !
  !---------!

  do p=nC+1,nBas-nV
    H(p,p) = eHF(p)
  end do

  !-------------!
  ! Block V2h1p !
  !-------------!

  do p=nC+1,nBas-nV

    klc = 0
    do k=nC+1,nO
      do l=nC+1,nO
        do c=nO+1,nBas-nV
          klc = klc + 1

          H(p       ,nBas+klc) = sqrt(2d0)*ERI(p,c,k,l)
          H(nBas+klc,p       ) = sqrt(2d0)*ERI(p,c,k,l)

        end do
      end do
    end do

  end do

  !-------------!
  ! Block V2p1h !
  !-------------!

  do p=nC+1,nBas-nV

    kcd = 0
    do k=nC+1,nO
      do c=nO+1,nBas-nV
        do d=nO+1,nBas-nV
          kcd = kcd + 1

          H(p       ,nBas+kcd) = sqrt(2d0)*ERI(p,k,d,c)
          H(nBas+kcd,p       ) = sqrt(2d0)*ERI(p,k,d,c)

        end do
      end do
    end do

  end do

  !-------------!
  ! Block C2h1p !
  !-------------!

  ija = 0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nV
        ija = ija + 1

        klc = 0
        do k=nC+1,nO
          do l=nC+1,nO
            do c=nO+1,nBas-nV
              klc = klc + 1

              H(nBas+ija,nBas+klc) = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
                                   - ERI(j,c,a,l))**Kronecker_delta(i,k)

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
    do a=nO+1,nBas-nV
      do b=nO+1,nBas-nV
        iab = iab + 1

        kcd = 0
        do k=nC+1,nO
          do c=nO+1,nBas-nV
            do d=nO+1,nBas-nV
              kcd = kcd + 1

              H(nBas+iab,nBas+kcd) = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                                   - ERI(j,c,a,l))**Kronecker_delta(b,d)

            end do
          end do
        end do

      end do
    end do
  end do

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

  call diagonalize_matrix(nH,H,eGW)

!--------------!
! Dump results !
!--------------!

  write(*,*) '---------------------------------------'
  write(*,*) ' GW supermatrix quasiparticle energies '
  write(*,*) '---------------------------------------'
  write(*,*) 
  call matout(nH,1,eGW)
  write(*,*) 

end subroutine ufGW
