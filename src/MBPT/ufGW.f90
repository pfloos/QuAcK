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
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: klc,kcd,ija,iab

  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)

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

  n2h1p = nO*nO*nS
  n2p1h = nV*nV*nO
  nH = nBas + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))

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

  do p=nC+1,nBas-nR
    H(p,p) = eHF(p)
  end do

  !-------------!
  ! Block V2h1p !
  !-------------!

  do p=nC+1,nBas-nR

    klc = 0
    do k=nC+1,nO
      do l=nC+1,nO
        do c=nO+1,nBas-nR
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

  do p=nC+1,nBas-nR

    kcd = 0
    do k=nC+1,nO
      do c=nO+1,nBas-nR
        do d=nO+1,nBas-nR
          kcd = kcd + 1

          H(p              ,nBas+n2h1P+kcd) = sqrt(2d0)*ERI(p,k,d,c)
          H(nBas+n2h1p+kcd,p              ) = sqrt(2d0)*ERI(p,k,d,c)

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
      do a=nO+1,nBas-nR
        ija = ija + 1

        klc = 0
        do k=nC+1,nO
          do l=nC+1,nO
            do c=nO+1,nBas-nR
              klc = klc + 1

              H(nBas+ija,nBas+klc) & 
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

              H(nBas+n2h1p+iab,nBas+n2h1p+kcd) &
                = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)

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

!-----------------!
! Compute weights !
!-----------------!

  Z(:) = 0d0
  do s=1,nH
    do p=nC+1,nBas-nR
      Z(s) = Z(s) + H(p,s)**2
    end do
  end do

!--------------!
! Dump results !
!--------------!

  write(*,*) '--------------------------------------------'
  write(*,*) ' GW supermatrix quasiparticle energies (eV) '
  write(*,*) '--------------------------------------------'
  write(*,*) 
  call matout(nH,1,HaToeV*eGW(:))
  write(*,*) 
  write(*,*) '-----------------------'
  write(*,*) ' Quasiparticle weights '
  write(*,*) '-----------------------'
  write(*,*) 
  call matout(nH,1,Z(:))

! write(*,*) 
! write(*,*) '-------------------------'
! write(*,*) ' GW supermatrix orbitals '
! write(*,*) '-------------------------'
! write(*,*) 
! call matout(nH,nH,H)
! write(*,*) 

end subroutine ufGW
