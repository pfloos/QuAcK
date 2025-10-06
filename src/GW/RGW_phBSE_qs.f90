subroutine RGW_phBSE_qs(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, & 
                        nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: ispin
  integer                       :: isp_W

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

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS),KA_sta(nS,nS),KB_sta(nS,nS), &
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

! Initialization

  EcBSE(:) = 0d0

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

                 call phRLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

               call RGW_phBSE_static_kernel_A_qs(eta,nOrb,nC,nO,nV,nR,nS,1d0,ERI,eGW,OmRPA,rho_RPA,KA_sta)
  if(.not.TDA) call RGW_phBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1

    ! Compute BSE excitation energies

                 call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phRLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','singlet',nS,OmBSE)
    call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2

    ! Compute BSE excitation energies

                 call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phRLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','triplet',nS,OmBSE)
    call phLR_transition_vectors(.false.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)


  end if

! Scale properly correlation energy if exchange is included in interaction kernel

  if(exchange_kernel) then

    EcBSE(1) = 0.5d0*EcBSE(1)
    EcBSE(2) = 1.5d0*EcBSE(2)

  end if

end subroutine RGW_phBSE_qs

subroutine RGW_phBSE_static_kernel_A_qs(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGW,Om,rho,KA)

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
             
            num = -rho(i,j,kc)*rho(a,b,kc)
             
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
