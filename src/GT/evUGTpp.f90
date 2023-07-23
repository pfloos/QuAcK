subroutine evUGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE, &
                   TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip,&
                   eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,ERI_aaaa, &
                   ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)

! Perform one-shot calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip 
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin) 
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond(nspin)
  double precision              :: Conv
  integer                       :: ispin,is
  integer                       :: iblock 
  integer                       :: nH_sc,nH_sf,nHaa,nHab,nHbb
  integer                       :: nP_sc,nP_sf,nPaa,nPab,nPbb
  double precision              :: EcRPA(nspin),Ecaa,Ecbb
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Om1ab(:),Om1aa(:),Om1bb(:)
  double precision,allocatable  :: X1ab(:,:),X1aa(:,:),X1bb(:,:)
  double precision,allocatable  :: Y1ab(:,:),Y1aa(:,:),Y1bb(:,:)
  double precision,allocatable  :: rho1ab(:,:,:),rho1aa(:,:,:),rho1bb(:,:,:)
  double precision,allocatable  :: Om2ab(:),Om2aa(:),Om2bb(:)
  double precision,allocatable  :: X2ab(:,:),X2aa(:,:),X2bb(:,:)
  double precision,allocatable  :: Y2ab(:,:),Y2aa(:,:),Y2bb(:,:)
  double precision,allocatable  :: rho2ab(:,:,:),rho2aa(:,:,:),rho2bb(:,:,:)
  double precision,allocatable  :: SigT(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: eGT(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: e_diis(:,:,:)
  

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Self-consistent evUGT calculation        |'
  write(*,*)'************************************************'
  write(*,*)

! Dimensions of the pp-URPA linear reponse matrices

  nPaa = nV(1)*(nV(1)-1)/2
  nPbb = nV(2)*(nV(2)-1)/2

  nHaa = nO(1)*(nO(1)-1)/2;
  nHbb = nO(2)*(nO(2)-1)/2;

  nPab = nV(1)*nV(2)
  nHab = nO(1)*nO(2)

  nP_sc = nPab
  nH_sc = nHab

  nP_sf = nPaa + nPbb
  nH_sf = nHaa + nHbb
  

! Memory allocation

  allocate(Om1ab(nPab),X1ab(nPab,nPab),Y1ab(nHab,nPab), & 
           Om2ab(nHab),X2ab(nPab,nHab),Y2ab(nHab,nHab), & 
           rho1ab(nBas,nBas,nPab),rho2ab(nBas,nBas,nHab), & 
           Om1aa(nPaa),X1aa(nPaa,nPaa),Y1aa(nHaa,nPaa), & 
           Om2aa(nHaa),X2aa(nPaa,nHaa),Y2aa(nHaa,nHaa), & 
           rho1aa(nBas,nBas,nPaa),rho2aa(nBas,nBas,nHaa), & 
           Om1bb(nPbb),X1bb(nPbb,nPbb),Y1bb(nHbb,nPbb), &
           Om2bb(nPbb),X2bb(nPbb,nPbb),Y2bb(nHbb,nPbb), &
           rho1bb(nBas,nBas,nPbb),rho2bb(nBas,nBas,nHbb), &
           SigT(nBas,nspin),Z(nBas,nspin), &
           eGT(nBas,nspin),eOld(nBas,nspin),error_diis(nBas,max_diis,nspin), &
           e_diis(nBas,max_diis,nspin))

!----------------------------------------------
! Compute the exchange part of the self-energy
!----------------------------------------------

!Initialization
  
  nSCF              = 0 
  n_diis            = 0
  Conv              = 1d0
  e_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0  
  eGT(:,:)          = eHF(:,:)
  eOld(:,:)         = eGT(:,:)
  Z(:,:)            = 1d0
  rcond(:)          = 0d0 

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

    ispin  = 1
    iblock = 3
! iblock = 1

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPab,nHaa,nHab,nHbb,nHab,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin)) 
  
!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

    ispin  = 2
    iblock = 4

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPaa,nHaa,nHab,nHbb,nHaa,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))  
  
!----------------------------------------------
! beta-beta block
!----------------------------------------------

    ispin  = 2
    iblock = 7

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPbb,nHaa,nHab,nHbb,nHbb,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1bb,X1bb,Y1bb,Om2bb,X2bb,Y2bb,EcRPA(ispin))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

    EcGM    = 0d0
    SigT(:,:) = 0d0
    Z(:,:)    = 0d0

!alpha-beta block
  
    iblock = 3

    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHab,nPab,ERI_aaaa,ERI_aabb,ERI_bbbb,X1ab,Y1ab, &
                                  rho1ab,X2ab,Y2ab,rho2ab)
!alpha-alpha block

    iblock = 4
  
    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHaa,nPaa,ERI_aaaa,ERI_aabb,ERI_bbbb,X1aa,Y1aa, &
                                  rho1aa,X2aa,Y2aa,rho2aa)

!beta-beta block 
  
    iblock = 7

    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHbb,nPbb,ERI_aaaa,ERI_aabb,ERI_bbbb,X1bb,Y1bb, &
                                  rho1bb,X2bb,Y2bb,rho2bb)

    call UGTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,nPab,nPbb,eGT,Om1aa,Om1ab,Om1bb,&
                                rho1aa,rho1ab,rho1bb,Om2aa,Om2ab,Om2bb,rho2aa,rho2ab,rho2bb,EcGM,SigT)

    call UGTpp_renormalization_factor(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,nPab,nPbb,eGT,Om1aa,Om1ab,&
                                      Om1bb,rho1aa,rho1ab,rho1bb,Om2aa,Om2ab,Om2bb,rho2aa,rho2ab,rho2bb,Z) 


    Z(:,:) = 1d0/(1d0 - Z(:,:))

!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------
  
    eGT(:,:) = eHF(:,:) + SigT(:,:) 

! Convergence criteria

    Conv = maxval(abs(eGT(:,:) - eOld(:,:)))

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

!alpha-beta block

    ispin  = 1
    iblock = 3 

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
               nPab,nHaa,nHab,nHbb,nHab,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1ab,X1ab,Y1ab, &
               Om2ab,X2ab,Y2ab,EcRPA(ispin))

!alpha-alpha block
 
    ispin  = 2
    iblock = 4 
  
    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
               nPaa,nHaa,nHab,nHbb,nHaa,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1aa,X1aa,Y1aa, &
               Om2aa,X2aa,Y2aa,EcRPA(ispin)) 
 
    Ecaa = EcRPA(2)

!beta-beta block
 
    iblock = 7
  
    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
               nPbb,nHaa,nHab,nHbb,nHbb,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1bb,X1bb,Y1bb, &
               Om2bb,X2bb,Y2bb,EcRPA(ispin))
  
    Ecbb = EcRPA(2) 
    EcRPA(2) = Ecaa + Ecbb
    EcRPA(1) = EcRPA(1) - EcRPA(2)
    EcRPA(2) = 3d0*EcRPA(2)

    call print_evUGT(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigT,Z,eGT,EcGM,EcRPA)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    do is=1,nspin
      call DIIS_extrapolation(rcond(ispin),nBas,nBas,n_diis, &
                              error_diis(:,1:n_diis,is), &
                              e_diis(:,1:n_diis,is),eGT(:,is)-eOld(:,is), &
                              eGT(:,is))
    end do

    ! Reset DIIS if required

    if(minval(rcond(:)) < 1d-15) n_diis = 0 
  
! Save quasiparticles energy for next cycle

    eOld(:,:) = eGT(:,:)

 ! Increment

    nSCF = nSCF + 1

  enddo

!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Free memory

    deallocate(Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab, &
               Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, &
               Om1bb,X1bb,Y1bb,Om2bb,X2bb,Y2bb,rho1bb,rho2bb)

end subroutine 
