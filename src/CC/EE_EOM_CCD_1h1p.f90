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
  double precision,allocatable  :: VL(:,:)
  double precision,allocatable  :: VR(:,:)
! double precision,allocatable  :: Leom(:,:,:)
! double precision,allocatable  :: Reom(:,:,:)

! integer                       :: nstate,m
! double precision              :: Ex,tmp

  integer,allocatable           :: order(:)

! double precision,allocatable  :: rdm1_oo(:,:)
! double precision,allocatable  :: rdm1_vv(:,:)

! double precision,allocatable  :: rdm2_oovv(:,:,:,:)
! double precision,allocatable  :: rdm2_ovvo(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************'
  write(*,*)'| EE-EOM-CCD (1h1p) |'
  write(*,*)'*********************'
  write(*,*)

! Size of the EOM Hamiltonian

  nS = (nO-nC)*(nV-nR)

! Memory allocation

  allocate(Foo(nO,nO),Fvv(nV,nV),Wovvo(nO,nV,nV,nO),H(nS,nS),Om(nS))
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

  allocate(VL(nS,nS),VR(nS,nS))

  if(nS > 0) then 

    call diagonalize_general_matrix_LR(nS,H,Om,VL,VR)

    do ia=1,nS
      order(ia) = ia
    end do

    call quick_sort(Om,order,nS)
    call set_order_LR(VL,VR,order,nS,nS)

    call print_excitation_energies('EE-EOM-CCD','spinorbital',nS,Om)

!   write(*,*) 'Right Eigenvectors'
!   call matout(nS,nS,VR)

!   call matout(nS,3,VR(:,1:3))

  end if

! allocate(Leom(nO,nV,nS),Reom(nO,nV,nS))

! do m=1,nS
!   ia = 0
!   do i=1,nO
!     do a=1,nV
!       ia = ia + 1
!       Leom(i,a,m) = VL(ia,m)
!       Reom(i,a,m) = VR(ia,m)
!     end do
!   end do
! end do

! deallocate(VL,VR)

!------------------------------------------------------------------------
! EOM section
!------------------------------------------------------------------------

! allocate(rdm1_oo(nO,nO),rdm1_vv(nV,nV))
! allocate(rdm2_oovv(nO,nO,nV,nV),rdm2_ovvo(nO,nV,nV,nO))

! nstate = 1

! tmp = 0d0
! do i=1,nO
!   do a=1,nV
!     tmp = tmp + Leom(i,a,nstate)*Reom(i,a,nstate)
!   end do
! end do
! print*,tmp

! rdm1_oo(:,:) = 0d0
! do i=1,nO
!   do j=1,nO
!     do c=1,nV

!       rdm1_oo(i,j) = rdm1_oo(i,j) - Reom(i,c,nstate)*Leom(j,c,nstate)

!     end do
!   end do
! end do

! rdm1_vv(:,:) = 0d0
! do a=1,nV
!   do b=1,nV
!     do k=1,nO

!       rdm1_vv(a,b) = rdm1_vv(a,b) + Reom(k,b,nstate)*Leom(k,a,nstate)

!     end do
!   end do
! end do

! rdm2_ovvo(:,:,:,:) = 0d0
! do i=1,nO
!   do a=1,nV
!     do b=1,nV
!       do j=1,nO
! 
!         rdm2_ovvo(i,a,b,j) = Reom(i,b,nstate)*Leom(j,a,nstate)

!       end do
!     end do
!   end do
! end do

! rdm2_oovv(:,:,:,:) = 0d0
! do i=1,nO
!   do j=1,nO
!     do a=1,nV
!       do b=1,nV

!         do k=1,nO
!           do c=1,nV
! 
!             rdm2_oovv(i,j,a,b) = rdm2_oovv(i,j,a,b) & 
!                                + Reom(j,b,nstate)*t(k,i,c,a)*Leom(k,c,nstate) &
!                                - Reom(i,b,nstate)*t(k,j,c,a)*Leom(k,c,nstate) &
!                                - Reom(j,a,nstate)*t(k,i,c,b)*Leom(k,c,nstate) &
!                                + Reom(i,a,nstate)*t(k,j,c,b)*Leom(k,c,nstate)

!           end do
!         end do

!       end do
!     end do
!   end do
! end do

! Ex = 0d0

! do i=1,nO
!   Ex = Ex + rdm1_oo(i,i)*eO(i)
! end do

! do a=1,nV
!   Ex = Ex + rdm1_vv(a,a)*eV(a)
! end do

! do i=1,nO
!   do a=1,nV
!     do b=1,nV
!       do j=1,nO
! 
!         Ex = Ex + rdm2_ovvo(i,a,b,j)*OVVO(i,a,b,j) + 0.25d0*rdm2_oovv(i,j,a,b)*OOVV(i,j,a,b)
!     
!       end do                  
!     end do                    
!   end do                      
! end do

! print*,'Ex = ',Ex
! print*,'Om = ',Om(nstate)

end subroutine 
