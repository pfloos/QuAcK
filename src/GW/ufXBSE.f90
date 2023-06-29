subroutine ufXBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF,OmRPA,sERI)

! Unfolded BSE+ equations

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
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: sERI(nBas,nBas,nS)

! Local variables

  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb,kc,iajb,kcld
  integer,parameter             :: maxH = 20
  double precision              :: eps1,eps2
  double precision              :: Ve,Vh,C2h2p

  integer                       :: n1h1p,n2h2p,nH
  double precision,external     :: Kronecker_delta

  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**********************************************'
  write(*,*)'|        Unfolded BSE+ calculation           |'
  write(*,*)'**********************************************'
  write(*,*)

! TDA for W

  write(*,*) 'Tamm-Dancoff approximation by default!'
  write(*,*)

! Dimension of the supermatrix

  n1h1p = nO*nV
  n2h2p = nO*nO*nV*nV
     nH = n1h1p + n2h2p

! Memory allocation

  allocate(H(nH,nH),X(nH,nH),Om(nH),Z(nH))

! Initialization

  H(:,:) = 0d0

!---------------------------!
!  Compute BSE+ supermatrix !
!---------------------------!
!                           !
!     |  A     Ve-Vh |      !
! H = |              |      ! 
!     | Ve-Vh  C2h2p |      ! 
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

        H(ia,jb) = (eHF(a) - eHF(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                 + 2d0*ERI(i,b,a,j) - ERI(i,b,j,a)

        do kc=1,nS

          do l=nC+1,nO
            eps1 = 1d0/(eHF(a) - eHF(l) + OmRPA(kc))
            eps2 = 1d0/(eHF(b) - eHF(l) + OmRPA(kc))
            H(ia,jb) = H(ia,jb) + Kronecker_delta(i,j)*sERI(a,l,kc)*sERI(b,l,kc)*(eps1+eps2)
          enddo

          do d=nO+1,nBas-nR
            eps1 = 1d0/(- eHF(i) + eHF(d) + OmRPA(kc))
            eps2 = 1d0/(- eHF(j) + eHF(d) + OmRPA(kc))
            H(ia,jb) = H(ia,jb) + Kronecker_delta(a,b)*sERI(i,d,kc)*sERI(j,d,kc)*(eps1+eps2)
          enddo

           eps1 = 1d0/(eHF(a) - eHF(i) + OmRPA(kc))
           eps2 = 1d0/(eHF(b) - eHF(j) + OmRPA(kc))
           H(ia,jb) = H(ia,jb) - 2d0*sERI(i,a,kc)*sERI(j,b,kc)*(eps1+eps2)

        end do

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

            Ve = sqrt(2d0)*Kronecker_delta(k,j)*ERI(b,a,c,i)
            Vh = sqrt(2d0)*Kronecker_delta(b,c)*ERI(a,k,i,j)

            H(n1h1p+iajb,kc        ) = Ve - Vh
            H(kc        ,n1h1p+iajb) = Ve - Vh
           
            end do
          end do

        end do
      end do
    end do
  end do

! iajb=0
! ia = 0
! do i=nC+1,nO
!   do a=nO+1,nBas-nR
!     ia = ia + 1
!     do j=nC+1,nO
!       do b=nO+1,nBas-nR
!       iajb = iajb + 1

!       kc = 0
!       do k=nC+1,nO
!         do c=nO+1,nBas-nR
!           kc = kc + 1

!           Ve = sqrt(2d0)*Kronecker_delta(k,j)*sERI(b,c,ia)
!           Vh = sqrt(2d0)*Kronecker_delta(b,c)*sERI(k,j,ia)

!           H(n1h1p+iajb,kc        ) = Ve - Vh
!           H(kc        ,n1h1p+iajb) = Ve - Vh
!          
!           end do
!         end do

!       end do
!     end do
!   end do
! end do

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

                  C2h2p = ((eHF(a) + eHF(b) - eHF(i) - eHF(j))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                        + 2d0*ERI(a,k,i,c))*Kronecker_delta(j,l)*Kronecker_delta(b,d)

                  H(n1h1p+iajb,n1h1p+kcld) = C2h2p

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

! iajb = 0
! ia = 0
! do i=nC+1,nO
!   do a=nO+1,nBas-nR
!     ia = ia + 1
!     do j=nC+1,nO
!       do b=nO+1,nBas-nR
!         iajb = iajb + 1

!         H(n1h1p+iajb,n1h1p+iajb) = Om(ia) + eHF(b) - eHF(j) 

!       end do
!     end do
!   end do
! end do

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

  X(:,:) = H(:,:)
  call diagonalize_matrix(nH,X,Om)

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
  write(*,*)'       BSE+ excitation energies (eV)       '
  write(*,*)'-------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','Omega (eV)','|','Z','|'
  write(*,*)'-------------------------------------------'

  do s=1,min(nH,maxH)
    if(Z(s) > 1d-7) &
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',Om(s)*HaToeV,'|',Z(s),'|'
  enddo

  write(*,*)'-------------------------------------------'
  write(*,*)

end subroutine 
