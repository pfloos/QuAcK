program QuAcK

  implicit none
  include 'parameters.h'

  logical                       :: doRHF,doUHF,doMOM 
  logical                       :: doMP2,doMP3,doMP2F12
  logical                       :: doCCD,doCCSD,doCCSDT
  logical                       :: doCIS,doTDHF,doADC
  logical                       :: doGF2,doGF3
  logical                       :: doG0W0,doevGW,doqsGW
  logical                       :: doMCMP2,doMinMCMP2
  logical                       :: doeNcusp
  integer                       :: nNuc,nBas,nBasCABS
  integer                       :: nEl(nspin),nC(nspin),nO(nspin),nV(nspin),nR(nspin),nS(nspin)
  double precision              :: ENuc,ERHF,EUHF,Norm
  double precision              :: EcMP2(3),EcMP3,EcMP2F12(3),EcMCMP2(3),Err_EcMCMP2(3),Var_EcMCMP2(3)

  double precision,allocatable  :: ZNuc(:),rNuc(:,:)
  double precision,allocatable  :: cHF(:,:,:),eHF(:,:),PHF(:,:,:)
  double precision,allocatable  :: eG0W0(:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:),KShell(:)
  double precision,allocatable  :: CenterShell(:,:),DShell(:,:),ExpShell(:,:)

  integer                       :: TrialType
  double precision,allocatable  :: cTrial(:),gradient(:),hessian(:,:)

  double precision,allocatable  :: S(:,:),T(:,:),V(:,:),Hc(:,:),X(:,:)
  double precision,allocatable  :: ERI_AO_basis(:,:,:,:),ERI_MO_basis(:,:,:,:)
  double precision,allocatable  :: F12(:,:,:,:),Yuk(:,:,:,:),FC(:,:,:,:,:,:)

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_MOM    ,end_MOM      ,t_MOM
  double precision              :: start_CCD    ,end_CCD      ,t_CCD
  double precision              :: start_CCSD   ,end_CCSD     ,t_CCSD
  double precision              :: start_CIS    ,end_CIS      ,t_CIS
  double precision              :: start_TDHF   ,end_TDHF     ,t_TDHF
  double precision              :: start_ADC    ,end_ADC      ,t_ADC
  double precision              :: start_GF2    ,end_GF2      ,t_GF2
  double precision              :: start_GF3    ,end_GF3      ,t_GF3
  double precision              :: start_G0W0   ,end_G0W0     ,t_G0W0
  double precision              :: start_evGW   ,end_evGW     ,t_evGW
  double precision              :: start_qsGW   ,end_qsGW     ,t_qsGW
  double precision              :: start_eNcusp ,end_eNcusp   ,t_eNcusp
  double precision              :: start_MP2    ,end_MP2      ,t_MP2
  double precision              :: start_MP3    ,end_MP3      ,t_MP3
  double precision              :: start_MP2F12 ,end_MP2F12   ,t_MP2F12
  double precision              :: start_MCMP2  ,end_MCMP2    ,t_MCMP2
  double precision              :: start_MinMCMP2,end_MinMCMP2,t_MinMCMP2

  integer                       :: maxSCF_HF,n_diis_HF
  double precision              :: thresh_HF
  logical                       :: DIIS_HF,guess_type,ortho_type

  integer                       :: maxSCF_CC,n_diis_CC
  double precision              :: thresh_CC
  logical                       :: DIIS_CC

  logical                       :: singlet_manifold,triplet_manifold

  integer                       :: maxSCF_GF,n_diis_GF,renormalization
  double precision              :: thresh_GF
  logical                       :: DIIS_GF

  integer                       :: maxSCF_GW,n_diis_GW
  double precision              :: thresh_GW
  logical                       :: DIIS_GW,COHSEX,SOSEX,BSE,TDA,G0W,GW0,linearize

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

! Which calculations do you want to do?

  call read_methods(doRHF,doUHF,doMOM,          &
                    doMP2,doMP3,doMP2F12,       &
                    doCCD,doCCSD,doCCSDT,       &
                    doCIS,doTDHF,doADC,         &
                    doGF2,doGF3,                &
                    doG0W0,doevGW,doqsGW,       &
                    doMCMP2)

! Read options for methods

  call read_options(maxSCF_HF,thresh_HF,DIIS_HF,n_diis_HF,guess_type,ortho_type,            &
                    maxSCF_CC,thresh_CC,DIIS_CC,n_diis_CC,                                  &
                    singlet_manifold,triplet_manifold,                                      &
                    maxSCF_GF,thresh_GF,DIIS_GF,n_diis_GF,renormalization,                  &
                    maxSCF_GW,thresh_GW,DIIS_GW,n_diis_GW,COHSEX,SOSEX,BSE,TDA,G0W,GW0,linearize, &
                    nMC,nEq,nWalk,dt,nPrint,iSeed,doDrift)

! Weird stuff

  doeNCusp   = .false.
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
  allocate(ZNuc(nNuc),rNuc(nNuc,3))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

  allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
           DShell(maxShell,maxK),ExpShell(maxShell,maxK))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nNuc,rNuc,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)
  nS(:) = nO(:)*nV(:)

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

  allocate(cHF(nBas,nBas,nspin),eHF(nBas,nspin),eG0W0(nBas),PHF(nBas,nBas,nspin), &
           S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas),     &
           ERI_AO_basis(nBas,nBas,nBas,nBas),ERI_MO_basis(nBas,nBas,nBas,nBas))

! Read integrals

  call read_integrals(nBas,S,T,V,Hc,ERI_AO_basis)  

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
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute RHF energy
!------------------------------------------------------------------------

  if(doUHF) then

    call cpu_time(start_HF)
    call UHF(maxSCF_HF,thresh_HF,n_diis_HF,guess_type,nBas,nO,S,T,V,Hc,ERI_AO_basis,X,ENuc,EUHF,eHF,cHF,PHF)
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,*)

  endif

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

  endif

!------------------------------------------------------------------------
! AO to MO integral transform for post-HF methods
!------------------------------------------------------------------------

  call AOtoMO_integral_transform(nBas,cHF,ERI_AO_basis,ERI_MO_basis)

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

  endif

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

  endif

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

  endif

!------------------------------------------------------------------------
! Perform CCD calculation
!------------------------------------------------------------------------

  if(doCCD) then

    call cpu_time(start_CCD)
    call CCD(maxSCF_CC,thresh_CC,n_diis_CC,nBas,nEl,ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CCD,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Perform CCSD or CCSD(T) calculation
!------------------------------------------------------------------------

  if(doCCSD) then

    call cpu_time(start_CCSD)
    call CCSD(maxSCF_CC,thresh_CC,n_diis_CC,doCCSDT,nBas,nEl,ERI_MO_basis,ENuc,ERHF,eHF)
    call cpu_time(end_CCSD)

    t_CCSD = end_CCSD - start_CCSD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCSD or CCSD(T)= ',t_CCSD,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CIS excitations
!------------------------------------------------------------------------

  if(doCIS) then

    call cpu_time(start_CIS)
    call CIS(singlet_manifold,triplet_manifold, & 
             nBas,nC,nO,nV,nR,nS,ERI_MO_basis,eHF)
    call cpu_time(end_CIS)

    t_CIS = end_CIS - start_CIS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS = ',t_CIS,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute TDHF excitations
!------------------------------------------------------------------------

  if(doTDHF) then

    call cpu_time(start_TDHF)
    call TDHF(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ERI_MO_basis,eHF)
    call cpu_time(end_TDHF)

    t_TDHF = end_TDHF - start_TDHF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for TDHF = ',t_TDHF,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute ADC excitations
!------------------------------------------------------------------------

  if(doADC) then

    call cpu_time(start_ADC)
    call ADC(singlet_manifold,triplet_manifold,maxSCF_GF,thresh_GF,n_diis_GF,nBas,nC,nO,nV,nR,eHF,ERI_MO_basis)
    call cpu_time(end_ADC)

    t_ADC = end_ADC - start_ADC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ADC = ',t_ADC,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute GF2 electronic binding energies
!------------------------------------------------------------------------

  if(doGF2) then

    call cpu_time(start_GF2)
!   call GF2(maxSCF_GF,thresh_GF,n_diis_GF,nBas,nC,nO,nV,nR,ERI_MO_basis,eHF)
    call GF2_diag(maxSCF_GF,thresh_GF,n_diis_GF,nBas,nC,nO,nV,nR,ERI_MO_basis,eHF)
    call cpu_time(end_GF2)

    t_GF2 = end_GF2 - start_GF2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF2,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute GF3 electronic binding energies
!------------------------------------------------------------------------

  if(doGF3) then

    call cpu_time(start_GF3)
    call GF3_diag(maxSCF_GF,thresh_GF,n_diis_GF,renormalization,nBas,nC,nO,nV,nR,ERI_MO_basis,eHF)
    call cpu_time(end_GF3)

    t_GF3 = end_GF3 - start_GF3
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF3 = ',t_GF3,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Perform G0W0 calculatiom
!------------------------------------------------------------------------

  eG0W0(:) = eHF(:,1)

  if(doG0W0) then
    
    call cpu_time(start_G0W0)
    call G0W0(COHSEX,SOSEX,BSE,TDA,singlet_manifold,triplet_manifold, & 
              nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,PHF,ERI_AO_basis,ERI_MO_basis,cHF,eHF,eG0W0)
    call cpu_time(end_G0W0)
  
    t_G0W0 = end_G0W0 - start_G0W0
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0W0 = ',t_G0W0,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Perform evGW calculation
!------------------------------------------------------------------------

  if(doevGW) then

    call cpu_time(start_evGW)
    call evGW(maxSCF_GW,thresh_GW,n_diis_GW, &
              COHSEX,SOSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold,linearize, &
              nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,ERI_AO_basis,ERI_MO_basis,PHF,cHF,eHF,eHF)
    call cpu_time(end_evGW)

    t_evGW = end_evGW - start_evGW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGW = ',t_evGW,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Perform qsGW calculation
!------------------------------------------------------------------------

  if(doqsGW) then 

    call cpu_time(start_qsGW)
    call qsGW(maxSCF_GW,thresh_GW,n_diis_GW, & 
              COHSEX,SOSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold, & 
              nBas,nC,nO,nV,nR,nS,ENuc,S,X,T,V,Hc,ERI_AO_basis,PHF,cHF,eHF)
    call cpu_time(end_qsGW)

    t_qsGW = end_qsGW - start_qsGW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGW = ',t_qsGW,' seconds'
    write(*,*)

  endif

!------------------------------------------------------------------------
! Compute e-N cusp dressing
!------------------------------------------------------------------------
  if(doeNcusp) then

    call cpu_time(start_eNcusp)
!   call eNcusp()
    call cpu_time(end_eNcusp)

    t_eNcusp = end_eNcusp - start_eNcusp
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for e-N cusp dressing  = ',t_eNcusp,' seconds'
    write(*,*)

  endif

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

  endif
!------------------------------------------------------------------------
! Compute MC-MP2 energy
!------------------------------------------------------------------------

  if(doMCMP2) then

    call cpu_time(start_MCMP2)
    call MCMP2(doDrift,nBas,nC,nO,nV,cHF,eHF,EcMP2,        &
               nMC,nEq,nWalk,dt,nPrint,                                  &
               nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
               Norm,EcMCMP2,Err_EcMCMP2,Var_EcMCMP2)
!   call MCMP2(.false.,doDrift,nBas,nEl,nC,nO,nV,cHF,eHF,EcMP2,        &
!              nMC,nEq,nWalk,dt,nPrint,                                  &
!              nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
!              TrialType,Norm,cTrial,gradient,hessian,                   &
!              EcMCMP2,Err_EcMCMP2,Var_EcMCMP2)
    call cpu_time(end_MCMP2)

    t_MCMP2 = end_MCMP2 - start_MCMP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MC-MP2 = ',t_MCMP2,' seconds'
    write(*,*)

  endif

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

  endif

!------------------------------------------------------------------------
! End of QuAcK
!------------------------------------------------------------------------
end program QuAcK
