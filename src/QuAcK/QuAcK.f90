program QuAcK

  implicit none
  include 'parameters.h'

  logical                       :: doRQuAcK,doUQuAcK,doGQuAcK
  logical                       :: doRHF,doUHF,doGHF,doROHF
  logical                       :: dostab,dosearch
  logical                       :: doMP2,doMP3
  logical                       :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical                       :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical                       :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical                       :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical                       :: doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3
  logical                       :: doG0W0,doevGW,doqsGW,doufG0W0,doufGW,doSRGqsGW
  logical                       :: doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh

  integer                       :: nNuc,nBas
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
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: dipole_int_AO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)

  double precision              :: start_QuAcK,end_QuAcK,t_QuAcK
  double precision              :: start_int  ,end_int  ,t_int

  integer                       :: maxSCF_HF,max_diis_HF
  double precision              :: thresh_HF,level_shift,mix
  integer                       :: guess_type

  logical                       :: reg_MP

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

  call read_methods(doRHF,doUHF,doGHF,doROHF,              &
                    doMP2,doMP3,                           &
                    doCCD,dopCCD,doDCD,doCCSD,doCCSDT,     &
                    dodrCCD,dorCCD,docrCCD,dolCCD,         &
                    doCIS,doCIS_D,doCID,doCISD,doFCI,      & 
                    dophRPA,dophRPAx,docrRPA,doppRPA,      &
                    doG0F2,doevGF2,doqsGF2,doufG0F02,      & 
                    doG0F3,doevGF3,                        &
                    doG0W0,doevGW,doqsGW,doSRGqsGW,        &
                    doufG0W0,doufGW,                       &
                    doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp, &
                    doG0T0eh,doevGTeh,doqsGTeh,            &
                    doRtest,doUtest,doGtest)

!--------------------------!
! Read options for methods !
!--------------------------!

  call read_options(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,dostab,dosearch, &
                    reg_MP,                                                                     &
                    maxSCF_CC,thresh_CC,max_diis_CC,                                            &
                    TDA,spin_conserved,spin_flip,                                               &
                    maxSCF_GF,thresh_GF,max_diis_GF,lin_GF,eta_GF,renorm_GF,reg_GF,             &
                    maxSCF_GW,thresh_GW,max_diis_GW,lin_GW,eta_GW,reg_GW,TDA_W,                 &  
                    maxSCF_GT,thresh_GT,max_diis_GT,lin_GT,eta_GT,reg_GT,TDA_T,                 & 
                    doACFDT,exchange_kernel,doXBS,                                              &
                    dophBSE,dophBSE2,doppBSE,dBSE,dTDA)

!------------------------------------------------!
! Read input information                         !
!------------------------------------------------!
! nC   = number of core orbitals                 !
! nO   = number of occupied orbitals             !
! nV   = number of virtual orbitals (see below)  !
! nR   = number of Rydberg orbitals              !
! nBas = number of basis functions (see below)   !
!------------------------------------------------!

  call read_molecule(nNuc,nO,nC,nR)
  allocate(ZNuc(nNuc),rNuc(nNuc,ncart))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

!---------------------------------------!
! Read basis set information from PySCF !
!---------------------------------------!

  call read_basis_pyscf(nBas,nO,nV)

!--------------------------------------!
! Read one- and two-electron integrals !
!--------------------------------------!

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), & 
           ERI_AO(nBas,nBas,nBas,nBas),dipole_int_AO(nBas,nBas,ncart))

! Read integrals

  call wall_time(start_int)

  call read_integrals(nBas,S,T,V,Hc,ERI_AO)
  call read_dipole_integrals(nBas,dipole_int_AO)

  call wall_time(end_int)

  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for reading integrals = ',t_int,' seconds'
  write(*,*)

! Compute orthogonalization matrix

  call orthogonalization_matrix(nBas,S,X)

!---------------------!
! Choose QuAcK branch !
!---------------------!

  doRQuAcK = .false.
  if(doRHF .or. doROHF) doRQuAcK = .true.

  doUQuAcK = .false.
  if(doUHF) doUQuAcK = .true.

  doGQuAcK = .false.
  if(doGHF) doGQuAcK = .true.

!-----------------!
! Initialize Test !
!-----------------!

  dotest = doRtest .or. doUtest .or. doGtest

  if(dotest) call init_test(doRtest,doUtest,doGtest)

!-------------------------!
! Restricted QuAcK branch !
!-------------------------!

  if(doRQuAcK) &
    call RQuAcK(doRtest,doRHF,doROHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,              &
                dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA, &
                doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,doSRGqsGW,  &
                doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                                & 
                nNuc,nBas,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                            &
                S,T,V,Hc,X,dipole_int_AO,ERI_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                     &
                guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,              &
                maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,  &
                TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,           &
                dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

!---------------------------!
! Unrestricted QuAcK branch !
!---------------------------!

  if(doUQuAcK) &
    call UQuAcK(doUtest,doUHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,                     &
                dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA, &
                doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,doSRGqsGW,  &
                doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                                & 
                nNuc,nBas,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                            &
                S,T,V,Hc,X,dipole_int_AO,ERI_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                     &
                guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,              &
                maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,  &
                TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,           &
                dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

!--------------------------!
! Generalized QuAcK branch !
!--------------------------!

  if(doGQuAcK) & 
    call GQuAcK(doGtest,doGHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,              &
                dodrCCD,dorCCD,docrCCD,dolCCD,dophRPA,dophRPAx,docrRPA,doppRPA,                           &  
                doG0W0,doevGW,doqsGW,doG0F2,doevGF2,doqsGF2,                                              &
                nNuc,nBas,sum(nC),sum(nO),sum(nV),sum(nR),ENuc,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,ERI_AO, &
                maxSCF_HF,max_diis_HF,thresh_HF,level_shift,guess_type,mix,reg_MP,                        & 
                maxSCF_CC,max_diis_CC,thresh_CC,TDA,maxSCF_GF,max_diis_GF,thresh_GF,lin_GF,reg_GF,eta_GF, & 
                maxSCF_GW,max_diis_GW,thresh_GW,TDA_W,lin_GW,reg_GW,eta_GW,                               &
                dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

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
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for QuAcK = ',t_QuAcK,' seconds'
  write(*,*)

end program 
