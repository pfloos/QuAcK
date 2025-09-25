subroutine RQuAcK(working_dir,use_gpu,dotest,doRHF,doROHF,docRHF,                                                        &
                  dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,                                         &
                  dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA,       & 
                  doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,                  &
                  doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,doevParquet,doqsParquet,              & 
                  docG0W0,docG0F2,                                                                                       & 
                  doCAP,readFCIDUMP,                                                                                     & 
                  nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                                             &
                  S,T,V,Hc,CAP_AO,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                           &
                  guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,singlet,triplet,TDA,                             &
                  maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,        & 
                  TDA_W,lin_GW,reg_GW,eta_GW,doOO,mu,do_linDM_GW,                                                        &
                  maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,do_linDM_GT,                                & 
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS,                                      &
                  TDAeh,TDApp,max_diis_1b,max_diis_2b,max_it_1b,conv_1b,max_it_2b,conv_2b,lin_parquet,reg_1b,reg_2b,reg_PA)

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: use_gpu

  logical,intent(in)            :: dotest
  logical,intent(in)            :: readFCIDUMP

  logical,intent(in)            :: doRHF,doROHF,docRHF
  logical,intent(in)            :: dostab
  logical,intent(in)            :: dosearch
  logical,intent(in)            :: doaordm
  logical,intent(in)            :: doMP2,doMP3
  logical,intent(in)            :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical,intent(in)            :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical,intent(in)            :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical,intent(in)            :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical,intent(in)            :: doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3
  logical,intent(in)            :: doG0W0,doevGW,doqsGW,doufG0W0,doufGW
  logical,intent(in)            :: doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp
  logical,intent(in)            :: doG0T0eh,doevGTeh,doqsGTeh
  logical,intent(in)            :: docG0W0,docG0F2
  logical,intent(in)            :: doCAP
  logical,intent(in)            :: doevParquet,doqsParquet

  integer,intent(in)            :: nNuc,nBas,nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(inout):: ENuc

  double precision,intent(in)   :: ZNuc(nNuc),rNuc(nNuc,ncart)

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: CAP_AO(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

  logical,intent(in)            :: reg_MP

  integer,intent(in)            :: maxSCF_CC,max_diis_CC
  double precision,intent(in)   :: thresh_CC

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: TDA

  integer,intent(in)            :: maxSCF_GF,max_diis_GF,renorm_GF
  double precision,intent(in)   :: thresh_GF
  logical,intent(in)            :: lin_GF,reg_GF
  double precision,intent(in)   :: eta_GF

  integer,intent(in)            :: maxSCF_GW,max_diis_GW
  double precision,intent(in)   :: thresh_GW
  logical,intent(in)            :: TDA_W,lin_GW,reg_GW,doOO
  integer,intent(in)            :: mu
  double precision,intent(in)   :: eta_GW
  logical,intent(in)            :: do_linDM_GW

  integer,intent(in)            :: maxSCF_GT,max_diis_GT
  double precision,intent(in)   :: thresh_GT
  logical,intent(in)            :: TDA_T,lin_GT,reg_GT
  double precision,intent(in)   :: eta_GT
  logical,intent(in)            :: do_linDM_GT

  logical,intent(in)            :: dophBSE,dophBSE2,doppBSE,dBSE,dTDA
  logical,intent(in)            :: doACFDT,exchange_kernel,doXBS

  integer,intent(in)            :: max_it_1b,max_it_2b
  double precision,intent(in)   :: conv_1b,conv_2b
  integer,intent(in)            :: max_diis_1b,max_diis_2b
  logical,intent(in)            :: TDAeh,TDApp
  double precision,intent(in)   :: reg_1b,reg_2b
  logical,intent(in)            :: lin_parquet,reg_PA

! Local variables

  logical                       :: doMP,doCC,doCI,doRPA,doGF,doGW,doGT
  logical                       :: file_exists

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: start_MP     ,end_MP       ,t_MP
  double precision              :: start_CC     ,end_CC       ,t_CC
  double precision              :: start_CI     ,end_CI       ,t_CI
  double precision              :: start_RPA    ,end_RPA      ,t_RPA
  double precision              :: start_GF     ,end_GF       ,t_GF
  double precision              :: start_GW     ,end_GW       ,t_GW
  double precision              :: start_GT     ,end_GT       ,t_GT
  double precision              :: start_Parquet,end_Parquet  ,t_Parquet

  double precision              :: start_int, end_int, t_int
  double precision,allocatable  :: eHF(:)
  complex*16,allocatable        :: complex_eHF(:)
  double precision,allocatable  :: cHF(:,:)
  complex*16,allocatable        :: complex_cHF(:,:)
  double precision,allocatable  :: PHF(:,:)
  complex*16,allocatable        :: complex_PHF(:,:)
  double precision,allocatable  :: FHF(:,:)
  complex*16,allocatable        :: complex_FHF(:,:)
  double precision              :: ERHF
  double precision              :: Val
  complex*16                    :: complex_ERHF
  double precision,allocatable  :: CAP_MO(:,:)
  complex*16,allocatable        :: complex_CAP_MO(:,:)
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  complex*16,allocatable        :: complex_dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: ixyz
  integer                       :: nS
  complex*16,allocatable        :: complex_ERI_MO(:,:,:,:)
  double precision,allocatable  :: eGW(:)

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Restricted Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!
  if (docRHF) then
    allocate(complex_PHF(nBas,nBas))
    allocate(complex_eHF(nOrb))
    allocate(complex_cHF(nBas,nOrb))
    allocate(complex_FHF(nBas,nBas))
    allocate(complex_dipole_int_MO(nOrb,nOrb,ncart))
    allocate(dipole_int_MO(0,0,0))
    allocate(complex_ERI_MO(nOrb,nOrb,nOrb,nOrb))
    allocate(CAP_MO(0,0))
    if (doCAP) then 
      allocate(complex_CAP_MO(nOrb,nOrb))
    else
      allocate(complex_CAP_MO(0,0))
    end if
  else 
    allocate(PHF(nBas,nBas))
    allocate(eHF(nOrb))
    allocate(cHF(nBas,nOrb))
    allocate(FHF(nBas,nBas))
    allocate(dipole_int_MO(nOrb,nOrb,ncart))
    allocate(complex_dipole_int_MO(0,0,0))
    allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
    allocate(complex_CAP_MO(0,0))
    if (doCAP) then
        allocate(CAP_MO(nOrb,nOrb))
    else
        allocate(CAP_MO(0,0))
    end if
  end if

  allocate(eGW(nOrb))
  
  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)

! For the FCIDUMP read two-body integrals

  inquire(file='FCIDUMP', exist=file_exists)
  if(file_exists .and. readFCIDUMP) then
   write(*,*)
   write(*,*) 'Reading FCIDUMP two-body integrals'
   write(*,*)
   ERI_AO=0d0 
   open(unit=314, form='formatted', file='FCIDUMP', status='old')
   do
    read(314,*) Val,iorb,jorb,korb,lorb
    if(korb==lorb .and. lorb==0) then
     if(iorb==jorb .and. iorb==0) then
      exit
     endif
    else
     ERI_AO(iorb,jorb,korb,lorb)=Val
    endif
   enddo
  endif
  close(314)
  
  call wall_time(end_int)
  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

!---------------------!
! Hartree-Fock module !
!---------------------!

  if(doRHF) then

    call wall_time(start_HF)
    call RHF(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

  end if

  if(doROHF) then

    call wall_time(start_HF)
    call ROHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
              nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ROHF = ',t_HF,' seconds'
    write(*,*)

  end if

  if(docRHF) then
    call wall_time(start_HF)
    call cRHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,ENuc, &
             nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_ERHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for cRHF = ',t_HF,' seconds'
    write(*,*)

  end if

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

  call wall_time(start_AOtoMO)

  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)

  if (docRHF) then 
    
    ! Transform to complex MOs

    ! Read and transform dipole-related integrals
    do ixyz=1,ncart
      call complex_AOtoMO(nBas,nOrb,complex_cHF,dipole_int_AO(1,1,ixyz),complex_dipole_int_MO(1,1,ixyz))
    end do
    ! 4-index transform 
    call complex_AOtoMO_ERI_RHF(nBas,nOrb,complex_cHF,ERI_AO,complex_ERI_MO)
    ! Transform CAP integrals
    if (doCAP) then
            call complex_AOtoMO(nBas,nOrb,complex_cHF,CAP_AO,complex_CAP_MO)
            complex_CAP_MO = (0d0,1d0)*complex_CAP_MO
    end if
  else

    ! Transform to real MOs

    ! Read and transform dipole-related integrals
    do ixyz=1,ncart
      call AOtoMO(nBas,nOrb,cHF,dipole_int_AO(1,1,ixyz),dipole_int_MO(1,1,ixyz))
    end do
    ! 4-index transform 
    call AOtoMO_ERI_RHF(nBas,nOrb,cHF,ERI_AO,ERI_MO)
    ! Transform CAP integrals
    if (doCAP) call AOtoMO(nBas,nOrb,cHF,CAP_AO,CAP_MO)
  end if
  call wall_time(end_AOtoMO)

  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)

!-----------------------------------!
! Stability analysis of HF solution !
!-----------------------------------!

  nS = (nO - nC)*(nV - nR)

  if(dostab) then

    call wall_time(start_stab)
    call RHF_stability(nOrb,nC,nO,nV,nR,nS,eHF,ERI_MO)
    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

  if(dosearch) then

    call wall_time(start_stab)
    call RHF_search(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                    nBas,nOrb,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,X, & 
                    ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

!-----------------------!
! Moller-Plesset module !
!-----------------------!

  doMP = doMP2 .or. doMP3

  if(doMP) then

    call wall_time(start_MP)
    call RMP(dotest,doMP2,doMP3,reg_MP,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for MP = ',t_MP,' seconds'
    write(*,*)

  end if

!------------------------!
! Coupled-cluster module !
!------------------------!

  doCC = doCCD .or. dopCCD .or. doDCD .or. doCCSD .or. doCCSDT .or. &  
         dodrCCD .or. dorCCD .or. docrCCD .or. dolCCD

  if(doCC) then

    call wall_time(start_CC)
    call RCC(dotest,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,dodrCCD,dorCCD,docrCCD,dolCCD,  & 
             maxSCF_CC,thresh_CC,max_diis_CC,nBas,nOrb,nC,nO,nV,nR,Hc,ERI_AO,ERI_MO, & 
             ENuc,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CC = ',t_CC,' seconds'
    write(*,*)

  end if

!----------------------------------!
! Configuration interaction module !
!----------------------------------!

  doCI = doCIS .or. doCID .or. doCISD .or. doFCI

  if(doCI) then

    call wall_time(start_CI)
    call RCI(dotest,doCIS,doCIS_D,doCID,doCISD,doFCI,singlet,triplet,nOrb, &
             nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eHF,ERHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CI = ',t_CI,' seconds'
    write(*,*)

  end if

!-----------------------------------!
! Random-phase approximation module !
!-----------------------------------!

  doRPA = dophRPA .or. dophRPAx .or. docrRPA .or. doppRPA

  if(doRPA) then

    call wall_time(start_RPA)
    call RRPA(use_gpu,dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet, &
              nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!-------------------------!
! Green's function module !
!-------------------------!

doGF = doG0F2 .or. doevGF2 .or. doqsGF2 .or. doufG0F02 .or. doG0F3 .or. doevGF3 .or. docG0F2

  if(doGF .and. .not. docRHF) then
    call wall_time(start_GF)
    call RGF(dotest,doG0F2,doevGF2,doqsGF2,doufG0F02,doG0F3,doevGF3,renorm_GF,maxSCF_GF, &
             thresh_GF,max_diis_GF,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,lin_GF,         &
             eta_GF,reg_GF,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,complex_ERHF,                    &
             S,X,T,V,Hc,ERI_AO,ERI_MO,CAP_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!---------------------------------!
! complex Green's function module !
!---------------------------------!

  if(doGF .and. docRHF) then
    call wall_time(start_GF)
    call complex_RGF(dotest,docG0F2,doevGF2,doqsGF2,maxSCF_GF,                                   &
             thresh_GF,max_diis_GF,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,lin_GF,         &
             eta_GF,reg_GF,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,complex_ERHF,            &
             S,X,T,V,Hc,ERI_AO,complex_ERI_MO,dipole_int_AO,complex_dipole_int_MO,&
             complex_PHF,complex_cHF,complex_eHF,CAP_AO, complex_CAP_MO)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)
  end if

!-----------!
! GW module !
!-----------!

  doGW = doG0W0 .or. doevGW .or. doqsGW .or. doufG0W0 .or. doufGW .or. docG0W0
  if(doGW .and. .not. docRHF) then
    
    call wall_time(start_GW)
    call RGW(dotest,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,maxSCF_GW,thresh_GW,max_diis_GW,                     & 
             doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,singlet,triplet,      &
             lin_GW,eta_GW,reg_GW,doOO,mu,do_linDM_GW,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T, &
             V,Hc,ERI_AO,ERI_MO,CAP_MO,dipole_int_AO,dipole_int_MO,PHF,FHF,cHF,eHF,eGW)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GW = ',t_GW,' seconds'
    write(*,*)

  end if


!-------------------!
! complex GW module !
!-------------------!

  if(doGW .and. docRHF) then
    call wall_time(start_GW)
    call complex_RGW(dotest,docG0W0,doevGW,doqsGW,maxSCF_GW,thresh_GW,max_diis_GW,        & 
             doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,singlet,triplet, &
             lin_GW,eta_GW,reg_GW,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,               &
             V,Hc,ERI_AO,complex_ERI_MO,CAP_AO,complex_CAP_MO,dipole_int_AO,&
             complex_dipole_int_MO,complex_PHF,complex_cHF,complex_eHF)
    call wall_time(end_GW)

    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GW = ',t_GW,' seconds'
    write(*,*)
  end if
!-----------------!
! T-matrix module !
!-----------------!

  doGT = doG0T0pp .or. doevGTpp .or. doqsGTpp .or. doufG0T0pp .or. doG0T0eh .or. doevGTeh .or. doqsGTeh

  if(doGT) then
    
    call wall_time(start_GT)
    call RGT(dotest,doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,                 &
             maxSCF_GT,thresh_GT,max_diis_GT,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,  &
<<<<<<< HEAD
             TDA_T,TDA,dBSE,dTDA,singlet,triplet,lin_GT,eta_GT,reg_GT,do_linDM_GT,nNuc,ZNuc,rNuc,ENuc,&
             nBas,nOrb,nC,nO_,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,                   &
=======
             TDA_T,TDA,dBSE,dTDA,singlet,triplet,lin_GT,eta_GT,reg_GT,nNuc,ZNuc,rNuc,ENuc,            &
             nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,                   &
>>>>>>> ea2c44239f3f1523e94361a22f1c0697e23584d8
             dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------!
!     Parquet module     !
!------------------------!

  if(doevParquet) then
    call wall_time(start_Parquet)
    call R_evParquet(TDAeh,TDApp,max_diis_1b,max_diis_2b,lin_parquet,reg_1b,reg_2b,reg_PA,ENuc,max_it_1b,conv_1b,max_it_2b,conv_2b, & 
                  nOrb,nC,nO,nV,nR,nS,ERHF,eHF,ERI_MO)
    call wall_time(end_Parquet)
  
    t_Parquet = end_Parquet - start_Parquet
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Parquet module = ', t_Parquet, ' seconds'
    write(*,*)

  end if

  if(doqsParquet) then
    call wall_time(start_Parquet)
    call R_qsParquet(TDAeh,TDApp,max_diis_1b,max_diis_2b,reg_1b,reg_2b,reg_PA,ENuc,max_it_1b,conv_1b,max_it_2b,conv_2b, & 
         nBas,nOrb,nC,nO,nV,nR,nS,ERHF,PHF,cHF,eHF,S,X,T,V,Hc,ERI_AO,ERI_MO)           
    call wall_time(end_Parquet)
  
    t_Parquet = end_Parquet - start_Parquet
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Parquet module = ', t_Parquet, ' seconds'
    write(*,*)

  end if
  
! Memory deallocation

  if (allocated(eHF)) deallocate(eHF)
  if (allocated(cHF)) deallocate(cHF)
  if (allocated(PHF)) deallocate(PHF)
  if (allocated(FHF)) deallocate(FHF)
  if (allocated(dipole_int_MO)) deallocate(dipole_int_MO)
  if (allocated(ERI_MO)) deallocate(ERI_MO)
  if (allocated(complex_ERI_MO)) deallocate(complex_ERI_MO)
  if (allocated(ERI_AO)) deallocate(ERI_AO)
  if (allocated(CAP_MO)) deallocate(CAP_MO)
  if (allocated(complex_CAP_MO)) deallocate(complex_CAP_MO)
  if (allocated(complex_eHF)) deallocate(complex_eHF)
  if (allocated(complex_cHF)) deallocate(complex_cHF)
  if (allocated(complex_PHF)) deallocate(complex_PHF)
end subroutine
