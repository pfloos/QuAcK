subroutine EE_EOM_CCD_1h1p(nC,nO,nV,nR,eO,eV,OOVV,OVVO,t)

! EE-EOM-CCD calculation up to 1h1p

  implicit none

! Input variables

  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: t(nO,nO,nV,nV)
  
! Local variables

  integer                       :: a,b,c,d
  integer                       :: i,j,k,l
  integer                       :: ia,jb
  integer                       :: nS
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: Fvv(:,:)
  double precision,allocatable  :: Foo(:,:)
  double precision,allocatable  :: Wovvo(:,:,:,:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: Z(:,:)

  integer,allocatable           :: order(:)

! Hello world

  write(*,*)
  write(*,*)'*********************'
  write(*,*)'| EE-EOM-CCD (1h1p) |'
  write(*,*)'*********************'
  write(*,*)

! Size of the EOM Hamiltonian

  nS = (nO-nC)*(nV-nR)

! Memory allocation

  allocate(Foo(nO,nO),Fvv(nV,nV),Wovvo(nO,nV,nV,nO),H(nS,nS),Om(nS),Z(nS,nS))
  allocate(order(nS))


! Form one-body terms

  do a=1,nV-nR
    do b=1,nV-nR
 
      Fvv(a,b) = eV(a)*Kronecker_delta(a,b) 

      do i=1,nO-nC
        do j=1,nO-nC
          do c=1,nV-nR
    
!           Fvv(a,b) = Fvv(a,b) - 0.5d0*OOVV(i,j,b,c)*t(i,j,a,c)

          end do
        end do
      end do

    end do
  end do

  do i=1,nO-nC
    do j=1,nO-nC

      Foo(i,j) = eO(i)*Kronecker_delta(i,j)

      do k=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR

!           Foo(i,j) = Foo(i,j) + 0.5d0*OOVV(i,k,a,b)*t(j,k,a,b)

          end do
        end do
      end do

    end do
  end do

! Form two-body terms

  do i=1,nO-nC
    do b=1,nV-nR
      do a=1,nV-nR
        do j=1,nO-nC
 
          Wovvo(i,b,a,j) = OVVO(i,b,a,j)

          do k=1,nO-nC
            do c=1,nV-nR
    
              Wovvo(i,b,a,j) = Wovvo(i,b,a,j) + OOVV(i,k,a,c)*t(k,j,c,b)

            end do
          end do

        end do
      end do
    end do
  end do

! Form EOM Hamiltonian

  ia = 0
  do i=1,nO-nC
    do a=1,nV-nR
      ia = ia + 1

      jb = 0
      do j=1,nO-nC
        do b=1,nV-nR
          jb = jb + 1

          H(ia,jb) = Fvv(a,b)*Kronecker_delta(i,j) - Kronecker_delta(a,b)*Foo(i,j) + Wovvo(i,b,a,j)

        end do
      end do

    end do
  end do

! Diagonalize EOM Hamiltonian

  if(nS > 0) then 

    call diagonalize_general_matrix(nS,H,Om,Z)

    do ia=1,nS
      order(ia) = ia
    end do

    call quick_sort(Om,order,nS)
    call set_order(Z,order,nS,nS)

    call print_excitation('EE-EOM-CCD  ',3,nS,Om)

  end if

end subroutine EE_EOM_CCD_1h1p
