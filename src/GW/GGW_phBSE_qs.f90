subroutine GGW_phBSE_qs(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

  double precision,allocatable  :: W(:,:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS), &
           Aph(nS,nS),Bph(nS,nS),KA_sta(nS,nS),KB_sta(nS,nS), &
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  EcRPA = 0d0

                 call phGLR_A(dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phGLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call GGW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call GGW_phBSE_static_kernel_A_qs(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGW,OmRPA,rho_RPA,KA_sta)
  call GGW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta)


!-----!
! TDA !
!-----!

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

!---------------------------------!
! Compute BSE excitation energies !
!---------------------------------!

  EcBSE = 0d0

               call phGLR_A(dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
  if(.not.TDA) call phGLR_B(dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

               Aph(:,:) = Aph(:,:) + KA_sta(:,:)
  if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

  call phGLR(TDA,nS,Aph,Bph,EcBSE,OmBSE,XpY_BSE,XmY_BSE)

  call print_excitation_energies('phBSE@GW@GHF','spinorbital',nS,OmBSE)
  call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

end subroutine GGW_phBSE_qs

subroutine GGW_phBSE_static_kernel_A_qs(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGW,Om,rho,KA)

! Compute the BSE static kernel for the resonant block

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: dem
  double precision              :: num
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: KA(nS,nS)

  KA(:,:) = 0d0
  
! Compute static kernel

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          do kc=1,nS
             
            num = -rho(i,j,kc)*rho(a,b,kc)/2d0
             
            dem = eGW(a) - eGW(b) - Om(kc)
            KA(ia,jb) = KA(ia,jb) + num*dem/(dem**2 + eta**2)

            dem = eGW(b) - eGW(a) - Om(kc)
            KA(ia,jb) = KA(ia,jb) + num*dem/(dem**2 + eta**2)
             
            dem = eGW(j) - eGW(i) - Om(kc)
            KA(ia,jb) = KA(ia,jb) + num*dem/(dem**2 + eta**2)
             
            dem = eGW(i) - eGW(j) - Om(kc)
            KA(ia,jb) = KA(ia,jb) + num*dem/(dem**2 + eta**2)
            
          end do

        end do
      end do
    end do
  end do

end subroutine 
