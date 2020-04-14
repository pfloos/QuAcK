program QuAcK

  implicit none
  include 'parameters.h'

  logical                       :: doSph
  logical                       :: doRHF,doUHF,doMOM 
  logical                       :: doMP2,doMP3,doMP2F12
  logical                       :: doCCD,doCCSD,doCCSDT
  logical                       :: do_drCCD,do_rCCD,do_lCCD,do_pCCD
  logical                       :: doCIS,doRPA,doRPAx
  logical                       :: doppRPA,doADC
  logical                       :: doG0F2,doevGF2,doG0F3,doevGF3
  logical                       :: doG0W0,doevGW,doqsGW
  logical                       :: doG0T0,doevGT,doqsGT
  logical                       :: doMCMP2,doMinMCMP2
  logical                       :: doBas

  integer                       :: nNuc,nBas,nBasCABS
  integer                       :: nEl(nspin),nC(nspin),nO(nspin),nV(nspin),nR(nspin)
  integer                       :: nS(nspin)
  double precision              :: ENuc,ERHF,EUHF,Norm
  double precision              :: EcMP2(3),EcMP3,EcMP2F12(3),EcMCMP2(3),Err_EcMCMP2(3),Var_EcMCMP2(3)

  double precision,allocatable  :: ZNuc(:),rNuc(:,:)
  double precision,allocatable  :: cHF(:,:,:),eHF(:,:),PHF(:,:,:)

  double precision,allocatable  :: eG0W0(:)
  double precision,allocatable  :: eG0T0(:)

  logical                       :: doACFDT
  logical                       :: exchange_kernel
  logical                       :: doXBS

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)
  integer,allocatable           :: max_ang_mom(:)
  double precision,allocatable  :: min_exponent(:,:)
  double precision,allocatable  :: max_exponent(:)

  integer                       :: TrialType
  double precision,allocatable  :: cTrial(:),gradient(:),hessian(:,:)

  double precision,allocatable  :: S(:,:),T(:,:),V(:,:),Hc(:,:),H(:,:),X(:,:)
  double precision,allocatable  :: ERI_AO_basis(:,:,:,:)
  double precision,allocatable  :: ERI_MO_basis(:,:,:,:)
  double precision,allocatable  :: F12(:,:,:,:),Yuk(:,:,:,:),FC(:,:,:,:,:,:)

  double precision              :: start_QuAcK  ,end_QuAcK    ,t_QuAcK
  double precision              :: start_int    ,end_int      ,t_int  
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_MOM    ,end_MOM      ,t_MOM
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: start_CCD    ,end_CCD      ,t_CCD
  double precision              :: start_CCSD   ,end_CCSD     ,t_CCSD
  double precision              :: start_CIS    ,end_CIS      ,t_CIS
  double precision              :: start_RPA    ,end_RPA      ,t_RPA 
  double precision              :: start_RPAx   ,end_RPAx     ,t_RPAx
  double precision              :: start_ppRPA  ,end_ppRPA    ,t_ppRPA
  double precision              :: start_ADC    ,end_ADC      ,t_ADC
  double precision              :: start_GF2    ,end_GF2      ,t_GF2
  double precision              :: start_GF3    ,end_GF3      ,t_GF3
  double precision              :: start_G0W0   ,end_G0W0     ,t_G0W0
  double precision              :: start_evGW   ,end_evGW     ,t_evGW
  double precision              :: start_qsGW   ,end_qsGW     ,t_qsGW
  double precision              :: start_G0T0   ,end_G0T0     ,t_G0T0
  double precision              :: start_evGT   ,end_evGT     ,t_evGT
  double precision              :: start_qsGT   ,end_qsGT     ,t_qsGT
  double precision              :: start_MP2    ,end_MP2      ,t_MP2
  double precision              :: start_MP3    ,end_MP3      ,t_MP3
  double precision              :: start_MP2F12 ,end_MP2F12   ,t_MP2F12
  double precision              :: start_MCMP2  ,end_MCMP2    ,t_MCMP2
  double precision              :: start_MinMCMP2,end_MinMCMP2,t_MinMCMP2
  double precision              :: start_Bas    ,end_Bas      ,t_Bas

  integer                       :: maxSCF_HF,n_diis_HF
  double precision              :: thresh_HF
  logical                       :: DIIS_HF,guess_type,ortho_type

  integer                       :: maxSCF_CC,n_diis_CC
  double precision              :: thresh_CC
  logical                       :: DIIS_CC

  logical                       :: singlet_manifold
  logical                       :: triplet_manifold

  integer                       :: maxSCF_GF,n_diis_GF,renormGF
  double precision              :: thresh_GF
  logical                       :: DIIS_GF,linGF

  integer                       :: maxSCF_GW,n_diis_GW
  double precision              :: thresh_GW
  logical                       :: DIIS_GW,COHSEX,SOSEX,BSE,TDA,G0W,GW0,linGW
  double precision              :: eta

  integer                       :: nMC,nEq,nWalk,nPrint,iSeed
  double precision              :: dt
  logical                       :: doDrift

! Hello World

  write(*,*)
  write(*,*) '******************************************************************************************'
  write(*,*) '*            QuAcK                       QuAcK                         QuAcK             *'
  write(*,*) '*   __        __        __       __        __        __       __        __        __     *'
  write(*,*) '* <(o )___  <(o )___  <(o )___ <(o )___  <(o )___  <(o )___ <(o )___  <(o )___  <(o )___ *'
  write(*,*) '* ( ._> /   ( ._> /   ( ._> /  ( ._> /   ( ._> /   ( ._> /  ( ._> /   ( ._> /   ( ._> /  *'
  write(*,*) '*|--------------------------------------------------------------------------------------|*'
  write(*,*) '******************************************************************************************'
  write(*,*)

! Spherium calculation?
  
  doSph = .false.

  call cpu_time(start_QuAcK)

! Which calculations do you want to do?

  call read_methods(doRHF,doUHF,doMOM,                &
                    doMP2,doMP3,doMP2F12,             &
                    doCCD,doCCSD,doCCSDT,             &
                    do_drCCD,do_rCCD,do_lCCD,do_pCCD, &
                    doCIS,doRPA,doRPAx,               &    
                    doppRPA,doADC,                    &
                    doG0F2,doevGF2,doG0F3,doevGF3,    &
                    doG0W0,doevGW,doqsGW,             &
                    doG0T0,doevGT,doqsGT,             &
                    doMCMP2)

! Read options for methods

  call read_options(maxSCF_HF,thresh_HF,DIIS_HF,n_diis_HF,guess_type,ortho_type, &
                    maxSCF_CC,thresh_CC,DIIS_CC,n_diis_CC,                       &
                    singlet_manifold,triplet_manifold,                           &
                    maxSCF_GF,thresh_GF,DIIS_GF,n_diis_GF,linGF,renormGF,        &
                    maxSCF_GW,thresh_GW,DIIS_GW,n_diis_GW,                       & 
                    COHSEX,SOSEX,BSE,TDA,G0W,GW0,linGW,eta,                      &  
                    doACFDT,exchange_kernel,doXBS,                               &
                    nMC,nEq,nWalk,dt,nPrint,iSeed,doDrift)

! Weird stuff

  doMinMCMP2 = .false.

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nC   = number of core orbitals
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nR   = number of Rydberg orbitals 
! nBas = number of basis functions (see below)
!      = nO + nV
! nS   = number of single excitation 
!      = nO*nV

  call read_molecule(nNuc,nEl(:),nO(:),nC(:),nR(:))
  allocate(ZNuc(nNuc),rNuc(nNuc,ncart))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

  allocate(CenterShell(maxShell,ncart),TotAngMomShell(maxShell),KShell(maxShell),DShell(maxShell,maxK), & 
           ExpShell(maxShell,maxK),max_ang_mom(nNuc),min_exponent(nNuc,maxL+1),max_exponent(nNuc))


!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nNuc,rNuc,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell, & 
                  max_ang_mom,min_exponent,max_exponent)
  nS(:) = (nO(:) - nC(:))*(nV(:) - nR(:))

!------------------------------------------------------------------------
! Read auxiliary basis set information
!------------------------------------------------------------------------

! call ReadAuxBasis(nNuc,rNuc,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell)

! Compute the number of basis functions

! call CalcNBasis(nShell,TotAngMomShell,nA)

! Number of virtual orbitals in complete space

! nBasCABS = nA - nBas

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(cHF(nBas,nBas,nspin),eHF(nBas,nspin),eG0W0(nBas),eG0T0(nBas),PHF(nBas,nBas,nspin),          &
           S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),H(nBas,nBas),X(nBas,nBas), &
           ERI_AO_basis(nBas,nBas,nBas,nBas),ERI_MO_basis(nBas,nBas,nBas,nBas))

! Read integrals

  call cpu_time(start_int)

  if(doSph) then

    call read_integrals_sph(nEl(:),nBas,S,T,V,Hc,ERI_AO_basis)  

  else

    call system('./GoQCaml')
    call read_integrals(nEl(:),nBas,S,T,V,Hc,ERI_AO_basis)

  end if

  call cpu_time(end_int)

    t_int = end_int - start_int
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for reading integrals = ',t_int,' seconds'
    write(*,*)

! Compute orthogonalization matrix

  call orthogonalization_matrix(ortho_type,nBas,S,X)

!------------------------------------------------------------------------
! Compute RHF energy
!------------------------------------------------------------------------

  if(doRHF) then

    call cpu_time(start_HF)
    call RHF(maxSCF_HF,thresh_HF,n_diis_HF,guess_type,nBas,nO,S,T,V,Hc,ERI_AO_basis,X,ENuc,ERHF,eHF,cHF,PHF)
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RHF = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute RHF energy
!------------------------------------------------------------------------

  if(doUHF) then

    call cpu_time(start_HF)
    call UHF(maxSCF_HF,thresh_HF,n_diis_HF,guess_type,nBas,nO,S,T,V,Hc,ERI_AO_basis,X,ENuc,EUHF,eHF,cHF,PHF)
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UHF = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Maximum overlap method
!------------------------------------------------------------------------

  if(doMOM) then

    call cpu_time(start_MOM)
    call MOM(maxSCF_HF,thresh_HF,n_diis_HF, &
             nBas,nO,S,T,V,Hc,ERI_AO_basis,X,ENuc,ERHF,cHF,eHF,PHF)
    call cpu_time(end_MOM)

    t_MOM = end_MOM - start_MOM
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MOM = ',t_MOM,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! AO to MO integral transform for post-HF methods
!------------------------------------------------------------------------

! Compute Hartree Hamiltonian in the MO basis

  call Hartree_matrix_MO_basis(nBas,cHF,PHF,Hc,ERI_AO_basis,H)

  call cpu_time(start_AOtoMO)

  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)


  if(doSph) then

    ERI_MO_basis = ERI_AO_basis
    print*,'!!! MO = AO !!!'
    deallocate(ERI_AO_basis)

  else

    call AOtoMO_integral_transform(nBas,cHF,ERI_AO_basis,ERI_MO_basis)

  end if

  call cpu_time(end_AOtoMO)

  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)

!------------------------------------------------------------------------
! Compute MP2 energy
!------------------------------------------------------------------------

  if(doMP2) then

    call cpu_time(start_MP2)
    call MP2(nBas,nC,nO,nV,nR,ERI_MO_basis,ENuc,ERHF,eHF,EcMP2)
    call cpu_time(end_MP2)

    t_MP2 = end_MP2 - start_MP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP2,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute MP3 energy
!------------------------------------------------------------------------

  if(doMP3) then

    call cpu_time(start_MP3)
    call MP3(nBas,nEl,ERI_MO_basis,eHF,ENuc,ERHF)
    call cpu_time(end_MP3)

    t_MP3 = end_MP3 - start_MP3
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP3 = ',t_MP3,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute MP2-F12 energy
!------------------------------------------------------------------------

  if(doMP2F12) then

    call cpu_time(start_MP2F12)

!   Memory allocation for one- and two-electron integrals

    allocate(F12(nBas,nBas,nBas,nBas),Yuk(nBas,nBas,nBas,nBas),FC(nBas,nBas,nBas,nBas,nBas,nBas))

!   Read integrals

    call read_F12_integrals(nBas,S,ERI_AO_basis,F12,Yuk,FC)
    call MP2F12(nBas,nC,nO,nV,ERI_AO_basis,F12,Yuk,FC,ERHF,eHF,cHF)
    call cpu_time(end_MP2F12)

    t_MP2F12 = end_MP2F12 - start_MP2F12
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2-F12 = ',t_MP2F12,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform CCD calculation
!------------------------------------------------------------------------

  if(doCCD) then

    call cpu_time(start_CCD)
    call CCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nC(1),nO(1),nV(1),nR(1), & 
             ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CCD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform CCSD or CCSD(T) calculation
!------------------------------------------------------------------------

  if(doCCSDT) doCCSD = .true.

  if(doCCSD) then

    call cpu_time(start_CCSD)
    call CCSD(maxSCF_CC,thresh_CC,n_diis_CC,doCCSDT,nBas,nC(1),nO(1),nV(1),nR(1), & 
              ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCSD)

    t_CCSD = end_CCSD - start_CCSD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCSD or CCSD(T)= ',t_CCSD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform direct ring CCD calculation
!------------------------------------------------------------------------

  if(do_drCCD) then

    call cpu_time(start_CCD)
    call drCCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nC(1),nO(1),nV(1),nR(1), &
               ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for direct ring CCD = ',t_CCD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ring CCD calculation
!------------------------------------------------------------------------

  if(do_rCCD) then

    call cpu_time(start_CCD)
    call rCCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nC(1),nO(1),nV(1),nR(1), &
              ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ring CCD = ',t_CCD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ladder CCD calculation
!------------------------------------------------------------------------

  if(do_lCCD) then

    call cpu_time(start_CCD)
    call lCCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nC(1),nO(1),nV(1),nR(1), &
              ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ladder CCD = ',t_CCD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform pair CCD calculation
!------------------------------------------------------------------------

  if(do_pCCD) then

    call cpu_time(start_CCD)
    call pCCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nC(1),nO(1),nV(1),nR(1), & 
              ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pair CCD = ',t_CCD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CIS excitations
!------------------------------------------------------------------------

  if(doCIS) then

    call cpu_time(start_CIS)
    call CIS(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ERI_MO_basis,eHF)
    call cpu_time(end_CIS)

    t_CIS = end_CIS - start_CIS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS = ',t_CIS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(doRPA) then

    call cpu_time(start_RPA)
    call RPA(doACFDT,exchange_kernel,singlet_manifold,triplet_manifold,eta, & 
             nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO_basis,eHF)
    call cpu_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute RPAx (RPA with exchange) excitations
!------------------------------------------------------------------------

  if(doRPAx) then

    call cpu_time(start_RPAx)
    call RPAx(doACFDT,exchange_kernel,singlet_manifold,triplet_manifold,eta, & 
              nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO_basis,eHF)
    call cpu_time(end_RPAx)

    t_RPAx = end_RPAx - start_RPAx
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPAx = ',t_RPAx,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute pp-RPA excitations
!------------------------------------------------------------------------

  if(doppRPA) then

    call cpu_time(start_ppRPA)
    call ppRPA(singlet_manifold,triplet_manifold, & 
               nBas,nC(1),nO(1),nV(1),nR(1),ENuc,ERHF,ERI_MO_basis,eHF)
    call cpu_time(end_ppRPA)

    t_ppRPA = end_ppRPA - start_ppRPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_ppRPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ADC excitations
!------------------------------------------------------------------------

  if(doADC) then

    call cpu_time(start_ADC)
    call ADC(singlet_manifold,triplet_manifold,maxSCF_GF,thresh_GF,n_diis_GF, & 
             nBas,nC(1),nO(1),nV(1),nR(1),eHF,ERI_MO_basis)
    call cpu_time(end_ADC)

    t_ADC = end_ADC - start_ADC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ADC = ',t_ADC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute G0F2 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F2) then

    call cpu_time(start_GF2)
    call G0F2(linGF,nBas,nC(1),nO(1),nV(1),nR(1),ERI_MO_basis,eHF)
    call cpu_time(end_GF2)

    t_GF2 = end_GF2 - start_GF2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF2,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF2 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF2) then

    call cpu_time(start_GF2)
    call evGF2(maxSCF_GF,thresh_GF,n_diis_GF,linGF,nBas,nC(1),nO(1),nV(1),nR(1),ERI_MO_basis,eHF)
    call cpu_time(end_GF2)

    t_GF2 = end_GF2 - start_GF2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF2,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute G0F3 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F3) then

    call cpu_time(start_GF3)
    call G0F3(renormGF,nBas,nC(1),nO(1),nV(1),nR(1),ERI_MO_basis,eHF)
    call cpu_time(end_GF3)

    t_GF3 = end_GF3 - start_GF3
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF3 = ',t_GF3,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF3 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF3) then

    call cpu_time(start_GF3)
    call evGF3(maxSCF_GF,thresh_GF,n_diis_GF,renormGF,nBas,nC(1),nO(1),nV(1),nR(1),ERI_MO_basis,eHF)
    call cpu_time(end_GF3)

    t_GF3 = end_GF3 - start_GF3
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF3 = ',t_GF3,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform G0W0 calculatiom
!------------------------------------------------------------------------

  eG0W0(:) = eHF(:,1)

  if(doG0W0) then
    
    call cpu_time(start_G0W0)
    call G0W0(doACFDT,exchange_kernel,doXBS,COHSEX,SOSEX,BSE,TDA, & 
              singlet_manifold,triplet_manifold,linGW,eta, & 
              nBas,nC(1),nO(1),nV(1),nR(1),nS(1),ENuc,ERHF,Hc,H,ERI_MO_basis,PHF,cHF,eHF,eG0W0)
    call cpu_time(end_G0W0)
  
    t_G0W0 = end_G0W0 - start_G0W0
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0W0 = ',t_G0W0,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGW calculation
!------------------------------------------------------------------------

  if(doevGW) then

    call cpu_time(start_evGW)
    call evGW(maxSCF_GW,thresh_GW,n_diis_GW,doACFDT,exchange_kernel,doXBS,COHSEX,SOSEX,BSE,TDA,G0W,GW0, &
              singlet_manifold,triplet_manifold,linGW,eta,                    &
              nBas,nC(1),nO(1),nV(1),nR(1),nS(1),ENuc,ERHF,Hc,H,ERI_MO_basis,PHF,cHF,eHF,eG0W0)
    call cpu_time(end_evGW)

    t_evGW = end_evGW - start_evGW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGW = ',t_evGW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGW calculation
!------------------------------------------------------------------------

  if(doqsGW) then 

    call cpu_time(start_qsGW)
    call qsGW(maxSCF_GW,thresh_GW,n_diis_GW,doACFDT,exchange_kernel,doXBS,COHSEX,SOSEX,BSE,TDA,G0W,GW0, & 
              singlet_manifold,triplet_manifold,eta, & 
              nBas,nC(1),nO(1),nV(1),nR(1),nS(1),ENuc,ERHF,S,X,T,V,Hc,ERI_AO_basis,ERI_MO_basis,PHF,cHF,eHF)
    call cpu_time(end_qsGW)

    t_qsGW = end_qsGW - start_qsGW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGW = ',t_qsGW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform G0T0 calculatiom
!------------------------------------------------------------------------

  eG0T0(:) = eHF(:,1)

  if(doG0T0) then
    
    call cpu_time(start_G0T0)
    call G0T0(doACFDT,exchange_kernel,doXBS,BSE,TDA,       &
              singlet_manifold,triplet_manifold,linGW,eta, &  
              nBas,nC(1),nO(1),nV(1),nR(1),nS(1),ENuc,ERHF,ERI_MO_basis,eHF)

    call cpu_time(end_G0T0)
  
    t_G0T0 = end_G0T0 - start_G0T0
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0T0 = ',t_G0T0,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGT calculatiom
!------------------------------------------------------------------------

  if(doevGT) then
    
    call cpu_time(start_evGT)
    call evGT(maxSCF_GW,thresh_GW,n_diis_GW,doACFDT,exchange_kernel,doXBS,BSE,TDA,singlet_manifold,triplet_manifold, &
                eta,nBas,nC(1),nO(1),nV(1),nR(1),nS(1),ENuc,ERHF,ERI_MO_basis,eHF,eG0T0)
    call cpu_time(end_evGT)
  
    t_evGT = end_evGT - start_evGT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGT = ',t_evGT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Information for Monte Carlo calculations
!------------------------------------------------------------------------

  if(doMCMP2 .or. doMinMCMP2) then

!   Print simulation details

    write(*,'(A32)') '----------------------'
    write(*,'(A32,1X,I16)')    'Number of Monte Carlo   steps',nMC
    write(*,'(A32,1X,I16)')    'Number of equilibration steps',nEq
    write(*,'(A32,1X,I16)')    'Number of walkers',nWalk
    write(*,'(A32,1X,F16.10)') 'Initial time step',dt
    write(*,'(A32,1X,I16)')    'Frequency of ouput',nPrint
    write(*,'(A32,1X,I16)')    'Seed for random number generator',iSeed
    write(*,'(A32)') '----------------------'
    write(*,*)

!   Initialize random number generator

    call initialize_random_generator(iSeed)

!------------------------------------------------------------------------
!   Type of weight function
!------------------------------------------------------------------------
!   TrialType = 0 => HF density
!   TrialType = 1 => Custom one-electron function
!------------------------------------------------------------------------

    TrialType = 0
    allocate(cTrial(nBas),gradient(nBas),hessian(nBas,nBas))

  end if
!------------------------------------------------------------------------
! Compute MC-MP2 energy
!------------------------------------------------------------------------

  if(doMCMP2) then

    call cpu_time(start_MCMP2)
    call MCMP2(doDrift,nBas,nC,nO,nV,cHF,eHF,EcMP2,        &
               nMC,nEq,nWalk,dt,nPrint,                                  &
               nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
               Norm,EcMCMP2,Err_EcMCMP2,Var_EcMCMP2)
    call cpu_time(end_MCMP2)

    t_MCMP2 = end_MCMP2 - start_MCMP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MC-MP2 = ',t_MCMP2,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Minimize MC-MP2 variance
!------------------------------------------------------------------------

  if(doMinMCMP2) then

    call cpu_time(start_MinMCMP2)
!   call MinMCMP2(nBas,nEl,nC,nO,nV,cHF,eHF,EcMP2,                              & 
!                 nMC,nEq,nWalk,dt,nPrint,                                  &
!                 nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, & 
!                 TrialType,Norm,cTrial,gradient,hessian)
    call cpu_time(end_MinMCMP2)
    
    t_MinMCMP2 = end_MinMCMP2 - start_MinMCMP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MC-MP2 variance minimization = ',t_MinMCMP2,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Basis set correction
!------------------------------------------------------------------------

  doBas = .false.

  if(doBas) then

    call cpu_time(start_Bas)
    call basis_correction(nBas,nO,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                          ERI_MO_basis,eHF,cHF,PHF,eG0W0)
    call cpu_time(end_Bas)

    t_Bas = end_Bas - start_Bas
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for basis set correction = ',t_Bas,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! End of QuAcK
!------------------------------------------------------------------------

  call cpu_time(end_QuAcK)

  t_QuAcK = end_QuAcK - start_QuAcK
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for QuAcK = ',t_QuAcK,' seconds'
  write(*,*)

end program QuAcK
