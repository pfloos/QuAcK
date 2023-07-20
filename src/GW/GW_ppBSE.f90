subroutine GW_ppBSE(TDA_W,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies at the pp level

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

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

  integer                       :: ispin
  integer                       :: isp_W

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: nOO
  integer                       :: nVV

  double precision,allocatable  :: A_sta(:,:)
  double precision,allocatable  :: B_sta(:,:)
  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)


  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: Om1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)

  double precision,allocatable  :: Om2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

  double precision,allocatable  :: WB(:,:)
  double precision,allocatable  :: WC(:,:)
  double precision,allocatable  :: WD(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

                 call phLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,A_sta)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

  call phLR(TDA_W,nS,A_sta,B_sta,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call GW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KA_sta)
  call GW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    write(*,*) '****************'
    write(*,*) '*** Singlets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 1
    EcBSE(ispin) = 0d0

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV), &
             Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
             WB(nVV,nOO),WC(nVV,nVV),WD(nOO,nOO))

    if(.not.TDA) call static_screening_WB_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WB)
    call static_screening_WC_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WC)
    call static_screening_WD_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WD)

    ! Compute BSE excitation energies

    call linear_response_pp_BSE(ispin,TDA,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,1d0,eGW,ERI,WB,WC,WD, &
                                Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))

!   call print_excitation('pp-BSE (N+2)',ispin,nVV,Om1)
!   call print_excitation('pp-BSE (N-2)',ispin,nOO,Om2)

    call print_transition_vectors_pp(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,WB,WC,WD)

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    write(*,*) '****************'
    write(*,*) '*** Triplets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 2
    EcBSE(ispin) = 0d0

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV), &
             Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
             WB(nVV,nOO),WC(nVV,nVV),WD(nOO,nOO))

    if(.not.TDA) call static_screening_WB_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WB)
    call static_screening_WC_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WC)
    call static_screening_WD_pp(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,WD)

    ! Compute BSE excitation energies

    call linear_response_pp_BSE(ispin,TDA,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,1d0,eGW,ERI,WB,WC,WD, &
                                Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))

!   call print_excitation('pp-BSE (N+2)',ispin,nVV,Om1)
!   call print_excitation('pp-BSE (N-2)',ispin,nOO,Om2)

    call print_transition_vectors_pp(.false.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,WB,WC,WD)

  end if

end subroutine 
