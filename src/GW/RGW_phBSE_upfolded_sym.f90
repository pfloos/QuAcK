subroutine RGW_phBSE_upfolded_sym(ispin,nBas,nOrb,nC,nO,nV,nR,nS,ERI,rho,Om,eHF)

! Upfolded phBSE@GW (TDA only!)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: m,ia,jb,iam,kcm
  integer,parameter             :: maxH = 20
  double precision              :: Jph,Kph,C2h2p

  integer                       :: n1h,n1p,n1h1p,n2h2p,nH
  double precision,external     :: Kronecker_delta

  integer,allocatable           :: order(:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************************'
  write(*,*)'* Symmetrized Upfolded phBSE@GW Calculation *'
  write(*,*)'*********************************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation by default!'
  write(*,*)

! Dimension of the supermatrix

  n1h = nO
  n1p = nV

  n1h1p = n1h*n1p
  n2h2p = n1h1p*n1h1p
     nH = n1h1p + n2h2p 

! Memory allocation

  allocate(order(nH),H(nH,nH),OmBSE(nH),Z(nH))

! Initialization

  H(:,:) = 0d0

!-------------------------------!
!  Compute BSE supermatrix      !
!-------------------------------!
!                               !
!     |     A     Jph + Kph |   ! 
! H = |                     |   ! 
!     | Jph + Kph   C2h2p   |   ! 
!                               !
!-------------------------------!

  !----------------------!
  ! Block A for singlets !
  !----------------------!

  if(ispin == 1) then

    ia = 0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        ia = ia + 1
       
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nOrb-nR
          jb = jb + 1
       
            H(ia,jb) = (eHF(a) - eHF(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                     + 2d0*ERI(i,b,a,j) - ERI(i,b,j,a)
       
            do k=nC+1,nO
              do m=1,nS

                H(ia,jb) = H(ia,jb) & 
                         + Kronecker_delta(i,j)*rho(a,k,m)*rho(b,k,m)/(eHF(a) - eHF(k) + Om(m)) &
                         + Kronecker_delta(i,j)*rho(a,k,m)*rho(b,k,m)/(eHF(b) - eHF(k) + Om(m))

              end do
            end do
       
            do c=nO+1,nOrb-nR
              do m=1,nS

                H(ia,jb) = H(ia,jb) & 
                         - Kronecker_delta(a,b)*rho(i,c,m)*rho(j,c,m)/(eHF(i) - eHF(c) - Om(m)) &
                         - Kronecker_delta(a,b)*rho(i,c,m)*rho(j,c,m)/(eHF(j) - eHF(c) - Om(m))

              end do
            end do
       
          end do
        end do
 
      end do
    end do

  end if

  !----------------------!
  ! Block A for triplets !
  !----------------------!

  if(ispin == 2) then

    ia = 0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        ia = ia + 1
 
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nOrb-nR
            jb = jb + 1
       
            H(ia,jb) = (eHF(a) - eHF(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                     - ERI(i,b,j,a)

            do k=nC+1,nO
              do m=1,nS

                H(ia,jb) = H(ia,jb) & 
                         + Kronecker_delta(i,j)*rho(a,k,m)*rho(b,k,m)/(eHF(a) - eHF(k) + Om(m)) &
                         + Kronecker_delta(i,j)*rho(a,k,m)*rho(b,k,m)/(eHF(b) - eHF(k) + Om(m))

              end do
            end do
       
            do c=nO+1,nOrb-nR
              do m=1,nS

                H(ia,jb) = H(ia,jb) & 
                         - Kronecker_delta(a,b)*rho(i,c,m)*rho(j,c,m)/(eHF(i) - eHF(c) - Om(m)) &
                         - Kronecker_delta(a,b)*rho(i,c,m)*rho(j,c,m)/(eHF(j) - eHF(c) - Om(m))

              end do
            end do
       
          end do
        end do
       
      end do
    end do

  end if

  !------------------!
  ! Blocks Jph & Kph !
  !------------------!

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ia = ia + 1

      kcm = 0
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          do m=1,nS
            kcm = kcm + 1

            ! Jph + Kph

            Jph = - sqrt(2d0)*Kronecker_delta(a,c)*rho(i,k,m)
            Kph = + sqrt(2d0)*Kronecker_delta(i,k)*rho(a,c,m)
            H(ia       ,n1h1p+kcm) = Jph + Kph
            H(n1h1p+kcm,ia)        = Jph + Kph

          end do
        end do
      end do

    end do
  end do

  !-------------!
  ! Block 2h2p !
  !-------------!

  iam = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do m=1,nS
        iam = iam + 1

        C2h2p = Om(m) + eHF(a) - eHF(i) 
        H(n1h1p+iam,n1h1p+iam) = C2h2p

      end do
    end do
  end do

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

  call diagonalize_matrix(nH,H,OmBSE)

!-----------------!
! Compute weights !
!-----------------!

  Z(:) = 0d0
  do s=1,nH
    do ia=1,n1h1p
      Z(s) = Z(s) + H(ia,s)**2
    end do
  end do

!--------------!
! Dump results !
!--------------!

  write(*,*)'-------------------------------------------'
  write(*,*)'  Upfolded phBSE excitation energies (eV)  '
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
