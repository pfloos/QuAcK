program QuAcK

  implicit none
  include 'parameters.h'

  logical                       :: doRQuAcK,doUQuAcK,doGQuAcK,doBQuAcK
  logical                       :: doRHF,doUHF,doGHF,doROHF,doHFB,docRHF
  logical                       :: dostab,dosearch
  logical                       :: doMP2,doMP3
  logical                       :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical                       :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical                       :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical                       :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical                       :: doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3
  logical                       :: doG0W0,doevGW,doqsGW,doufG0W0,doufGW
  logical                       :: docG0W0
  logical                       :: doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh

  integer                       :: nNuc
  integer                       :: nBas
  integer                       :: nOrb
  integer                       :: nC(nspin)
  integer                       :: nO(nspin)
  integer                       :: nV(nspin)
  integer                       :: nR(nspin)
  double precision              :: ENuc

  double precision,allocatable  :: ZNuc(:),rNuc(:,:)

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:),X_tmp(:,:)
  double precision, allocatable :: CAP(:,:)
  double precision,allocatable  :: dipole_int_AO(:,:,:)
  double precision,allocatable  :: Uvec(:,:), Uval(:)

  double precision              :: start_QuAcK,end_QuAcK,t_QuAcK
  double precision              :: start_int  ,end_int  ,t_int

  integer                       :: maxSCF_HF,max_diis_HF
  double precision              :: thresh_HF,level_shift,mix
  integer                       :: guess_type

  double precision              :: eta_cap

  logical                       :: reg_MP

  logical                       :: switch_hpc
  logical                       :: use_gpu

  integer                       :: maxSCF_CC,max_diis_CC
  double precision              :: thresh_CC

  logical                       :: spin_conserved
  logical                       :: spin_flip
  logical                       :: TDA
  
  integer                       :: maxSCF_GF,max_diis_GF,renorm_GF
  double precision              :: thresh_GF
  logical                       :: lin_GF,reg_GF
  double precision              :: eta_GF

  integer                       :: maxSCF_GW,max_diis_GW
  double precision              :: thresh_GW
  logical                       :: TDA_W,lin_GW,reg_GW
  double precision              :: eta_GW

  integer                       :: maxSCF_GT,max_diis_GT
  double precision              :: thresh_GT
  logical                       :: TDA_T,lin_GT,reg_GT
  double precision              :: eta_GT

  logical                       :: dophBSE,dophBSE2,doppBSE,dBSE,dTDA
  logical                       :: doACFDT,exchange_kernel,doXBS

  logical                       :: dotest,doRtest,doUtest,doGtest

  logical                       :: chem_pot_hf
  logical                       :: restart_hfb
  double precision              :: temperature,sigma

  character(len=256)            :: working_dir

  ! Check if the right number of arguments is provided
  if(command_argument_count() < 1) then
    print *, "No working directory provided."
    stop
  else
    call get_command_argument(1, working_dir)
  endif

!-------------!
! Hello World !
!-------------!

  write(*,*)
  write(*,*) '******************************************************************************************'
  write(*,*) '*            QuAcK                       QuAcK                         QuAcK             *'
  write(*,*) '*   __        __        __       __        __        __       __        __        __     *'
  write(*,*) '* <(o )___  <(o )___  <(o )___ <(o )___  <(o )___  <(o )___ <(o )___  <(o )___  <(o )___ *'
  write(*,*) '* ( ._> /   ( ._> /   ( ._> /  ( ._> /   ( ._> /   ( ._> /  ( ._> /   ( ._> /   ( ._> /  *'
  write(*,*) '*|--------------------------------------------------------------------------------------|*'
  write(*,*) '******************************************************************************************'
  write(*,*)

!-----------------------!
! Starting QuAcK timing !
!-----------------------!

  call wall_time(start_QuAcK)

!------------------!
! Method selection !
!------------------!

  call read_methods(working_dir,                           &
                    doRHF,doUHF,doGHF,doROHF,doHFB,docRHF, &
                    doMP2,doMP3,                           &
                    doCCD,dopCCD,doDCD,doCCSD,doCCSDT,     &
                    dodrCCD,dorCCD,docrCCD,dolCCD,         &
                    doCIS,doCIS_D,doCID,doCISD,doFCI,      & 
                    dophRPA,dophRPAx,docrRPA,doppRPA,      &
                    doG0F2,doevGF2,doqsGF2,doufG0F02,      & 
                    doG0F3,doevGF3,                        &
                    doG0W0,doevGW,doqsGW,doufG0W0,doufGW,  &
                    doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp, &
                    doG0T0eh,doevGTeh,doqsGTeh,            &
                    docG0W0,                               &
                    doRtest,doUtest,doGtest)

!--------------------------!
! Read options for methods !
!--------------------------!

  call read_options(working_dir,                                                                &
                    maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,dostab,dosearch, &
                    reg_MP,                                                                     &
                    maxSCF_CC,thresh_CC,max_diis_CC,                                            &
                    TDA,spin_conserved,spin_flip,                                               &
                    maxSCF_GF,thresh_GF,max_diis_GF,lin_GF,eta_GF,renorm_GF,reg_GF,             &
                    maxSCF_GW,thresh_GW,max_diis_GW,lin_GW,eta_GW,reg_GW,TDA_W,                 &  
                    maxSCF_GT,thresh_GT,max_diis_GT,lin_GT,eta_GT,reg_GT,TDA_T,                 & 
                    doACFDT,exchange_kernel,doXBS,                                              &
                    dophBSE,dophBSE2,doppBSE,dBSE,dTDA,                                         &
                    temperature,sigma,chem_pot_hf,restart_hfb)

!------------------!
! Hardware         !
!------------------!

  call read_hpc_flags(working_dir,switch_hpc,use_gpu)

!------------------------------------!
! Read input information             !
!------------------------------------!
! nC   = number of core orbitals     !
! nO   = number of occupied orbitals !
! nV   = number of virtual orbitals  !
! nR   = number of Rydberg orbitals  !
! nBas = number of basis functions   !
! nOrb = number of orbitals          !
!------------------------------------!

  call read_molecule(working_dir,nNuc,nO,nC,nR)
  allocate(ZNuc(nNuc),rNuc(nNuc,ncart))

! Read geometry

  call read_geometry(working_dir,nNuc,ZNuc,rNuc,ENuc)

!---------------------------------------!
! Read basis set information from PySCF !
!---------------------------------------!

  call read_basis_pyscf(working_dir,nBas,nO,nV)

!--------------------------------------!
! Read one- and two-electron integrals !
!--------------------------------------!

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas))
  allocate(T(nBas,nBas))
  allocate(V(nBas,nBas))
  allocate(Hc(nBas,nBas))
  allocate(dipole_int_AO(nBas,nBas,ncart))
  allocate(CAP(nBas,nBas))
! Read integrals

  call wall_time(start_int)

  call read_1e_integrals(working_dir,nBas,S,T,V,Hc)
  call read_eta_cap(working_dir,eta_cap)
  if (docRHF .or. docG0W0) call read_CAP_integrals(nBas,CAP) ! Add different cases if needed
  CAP(:,:) = -eta_cap*CAP(:,:)
  call read_dipole_integrals(working_dir,nBas,dipole_int_AO)
  call wall_time(end_int)

  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 1e-integrals = ',t_int,' seconds'
  write(*,*)

! Compute orthogonalization matrix

  allocate(X_tmp(nBas,nBas))
  call orthogonalization_matrix(nBas,nOrb,S,X_tmp)
  allocate(X(nBas,nOrb))
  X(1:nBas,1:nOrb) = X_tmp(1:nBas,1:nOrb)
  deallocate(X_tmp)

!---------------------!
! Choose QuAcK branch !
!---------------------!

  doRQuAcK = .false.
  if(doRHF .or. doROHF .or. docRHF) doRQuAcK = .true.

  doUQuAcK = .false.
  if(doUHF) doUQuAcK = .true.

  doGQuAcK = .false.
  if(doGHF) doGQuAcK = .true.

  doBQuAcK = .false.
  if(doHFB) doBQuAcK = .true.

!-----------------!
! Initialize Test !
!-----------------!

  dotest = doRtest .or. doUtest .or. doGtest

  if(dotest) call init_test(working_dir,doRtest,doUtest,doGtest)

!-------------------------!
! Restricted QuAcK branch !
!-------------------------!

  if(doRQuAcK) then

    if(switch_hpc) then
      call RQuAcK_hpc(working_dir,use_gpu,doRtest,doRHF,doROHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT, &
                      dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA,        &
                      doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,                   &
                      doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                                       &
                      nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                              &
                      S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                                   &
                      guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,                     &
                      maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,         &
                      TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,                  &
                      dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)
    else
      call RQuAcK(working_dir,use_gpu,doRtest,doRHF,doROHF,docRHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT, &
                  dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA,        &
                  doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,                   &
                  doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                                       &
                  docG0W0,                                                                                                &
                  nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                              &
                  S,T,V,Hc,CAP,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                       &
                  guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,                     &
                  maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,         &
                  TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,                  &
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)
    endif
  endif

!---------------------------!
! Unrestricted QuAcK branch !
!---------------------------!

  if(doUQuAcK) &
    call UQuAcK(working_dir,doUtest,doUHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,         &
                dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA, &
                doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,            &
                doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                                & 
                nNuc,nBas,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                            &
                S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                            &
                guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,              &
                maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,  &
                TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,           &
                dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

!--------------------------!
! Generalized QuAcK branch !
!--------------------------!
  if(doGQuAcK) & 
    call GQuAcK(working_dir,doGtest,doGHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,      &
                dodrCCD,dorCCD,docrCCD,dolCCD,dophRPA,dophRPAx,docrRPA,doppRPA,                               &
                doG0W0,doevGW,doqsGW,doG0F2,doevGF2,doqsGF2,doG0T0pp,doevGTpp,doqsGTpp,                       &
                nNuc,nBas,sum(nC),sum(nO),sum(nV),sum(nR),ENuc,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,            &
                maxSCF_HF,max_diis_HF,thresh_HF,level_shift,guess_type,mix,reg_MP,                            &
                maxSCF_CC,max_diis_CC,thresh_CC,TDA,maxSCF_GF,max_diis_GF,thresh_GF,lin_GF,reg_GF,eta_GF,     &
                maxSCF_GW,max_diis_GW,thresh_GW,TDA_W,lin_GW,reg_GW,eta_GW,                                   &
                maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,                                   &
                dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)


!--------------------------!
! Bogoliubov QuAcK branch !
!--------------------------!
  if(doBQuAcK) & 
    call BQuAcK(working_dir,dotest,doHFB,nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                           &
                S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,guess_type,mix,          &
                temperature,sigma,chem_pot_hf,restart_hfb)

!-----------!
! Stop Test !
!-----------!

  if(dotest) call stop_test(doRtest,doUtest,doGtest)

!--------------!
! Running Test !
!--------------!

  if(dotest) call run_test(doRtest,doUtest,doGtest)

!--------------!
! End of QuAcK !
!--------------!

  call wall_time(end_QuAcK)

  t_QuAcK = end_QuAcK - start_QuAcK
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for QuAcK = ',t_QuAcK,' seconds'
  write(*,*)

! Memory deallocation
  if (allocated(rNuc)) deallocate(rNuc)
  if (allocated(Znuc)) deallocate(Znuc)
  if (allocated(T)) deallocate(T)
  if (allocated(V)) deallocate(V)
  if (allocated(Hc)) deallocate(Hc)
  if (allocated(dipole_int_AO)) deallocate(dipole_int_AO)
  if (allocated(S)) deallocate(S)
end program 
