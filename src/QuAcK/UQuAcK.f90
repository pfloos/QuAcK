subroutine UQuAcK(working_dir,dotest,doUHF,docUHF,doMOM,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,    &
                  dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA,        &
                  doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,                                             &
                  doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,doevParquet,doqsParquet,               &
                  doCAP,                                                                                                  &
                  readFCIDUMP,nNuc,nBas,nC,nO,nV,nR,nCVS,FC,ENuc,ZNuc,rNuc,                                               &
                  S,T,V,Hc,CAP_AO,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                            &
                  mom_occupations,writeMOs,                                                                               &
                  guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,spin_conserved,spin_flip,TDA,                     &
                  maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,         &
                  TDA_W,lin_GW,reg_GW,eta_GW,maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,                  &
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest
  logical,intent(in)            :: readFCIDUMP

  logical,intent(in)            :: doUHF,docUHF,doMOM,doCAP,writeMOs
  logical,intent(in)            :: dostab
  logical,intent(in)            :: dosearch
  logical,intent(in)            :: doMP2,doMP3
  logical,intent(in)            :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical,intent(in)            :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical,intent(in)            :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical,intent(in)            :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical,intent(in)            :: doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3
  logical,intent(in)            :: doG0W0,doevGW,doqsGW
  logical,intent(in)            :: doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp
  logical,intent(in)            :: doG0T0eh,doevGTeh,doqsGTeh
  logical,intent(in)            :: doevParquet,doqsParquet

  integer,intent(in)            :: nNuc,nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  double precision,intent(inout):: ENuc
  integer,intent(in)            :: mom_occupations(maxval(nO),nspin)

  double precision,intent(in)   :: ZNuc(nNuc),rNuc(nNuc,ncart)

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: CAP_AO(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

  logical,intent(in)            :: reg_MP

  integer,intent(in)            :: maxSCF_CC,max_diis_CC
  double precision,intent(in)   :: thresh_CC

  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: TDA

  integer,intent(in)            :: maxSCF_GF,max_diis_GF,renorm_GF
  double precision,intent(in)   :: thresh_GF
  logical,intent(in)            :: lin_GF,reg_GF
  double precision,intent(in)   :: eta_GF

  integer,intent(in)            :: maxSCF_GW,max_diis_GW
  double precision,intent(in)   :: thresh_GW
  logical,intent(in)            :: TDA_W,lin_GW,reg_GW
  double precision,intent(in)   :: eta_GW

  integer,intent(in)            :: maxSCF_GT,max_diis_GT
  double precision,intent(in)   :: thresh_GT
  logical,intent(in)            :: TDA_T,lin_GT,reg_GT
  double precision,intent(in)   :: eta_GT

  logical,intent(in)            :: dophBSE,dophBSE2,doppBSE,dBSE,dTDA
  logical,intent(in)            :: doACFDT,exchange_kernel,doXBS

! Local variables

  logical                       :: file_exists
  logical                       :: doMP,doCC,doCI,doRPA,doGF,doGW,doGT
  logical                       :: CVS
  
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
  double precision,allocatable  :: cHF(:,:,:),eHF(:,:),PHF(:,:,:),FHF(:,:,:)
  complex*16,allocatable        :: complex_eHF(:,:),complex_cHF(:,:,:),complex_PHF(:,:,:),complex_FHF(:,:,:)
  double precision              :: Val
  double precision              :: EUHF
  complex*16                    :: complex_EUHF
  double precision,allocatable  :: CAP_MO(:,:,:)
  complex*16,allocatable        :: complex_CAP_MO(:,:,:)
  double precision,allocatable  :: dipole_int_aa(:,:,:),dipole_int_bb(:,:,:)
  complex*16,allocatable        :: complex_dipole_int_aa(:,:,:),complex_dipole_int_bb(:,:,:)
  double precision,allocatable  :: ERI_aaaa(:,:,:,:),ERI_aabb(:,:,:,:),ERI_bbbb(:,:,:,:)
  complex*16,allocatable        :: complex_ERI_aaaa(:,:,:,:)
  complex*16,allocatable        :: complex_ERI_aabb(:,:,:,:)
  complex*16,allocatable        :: complex_ERI_bbbb(:,:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: ixyz
  integer                       :: nS(nspin)

  write(*,*)
  write(*,*) '********************************'
  write(*,*) '* Unrestricted Branch of QuAcK *'
  write(*,*) '********************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------! 
  if(docUHF) then
    allocate(complex_PHF(nBas,nBas,nspin))
    allocate(complex_eHF(nBas,nspin))
    allocate(complex_cHF(nBas,nBas,nspin))
    allocate(complex_FHF(nBas,nBas,nspin))
    allocate(complex_dipole_int_aa(nBas,nBas,ncart))
    allocate(complex_dipole_int_bb(nBas,nBas,ncart))
    allocate(dipole_int_aa(0,0,0))
    allocate(dipole_int_bb(0,0,0))
    allocate(complex_ERI_aaaa(nBas,nBas,nBas,nBas))
    allocate(complex_ERI_aabb(nBas,nBas,nBas,nBas))
    allocate(complex_ERI_bbbb(nBas,nBas,nBas,nBas))
    allocate(CAP_MO(0,0,0))
    if (doCAP) then 
      allocate(complex_CAP_MO(nBas,nBas,nspin))
    else
      allocate(complex_CAP_MO(0,0,0))
    endif
  else 
    allocate(PHF(nBas,nBas,nspin),FHF(nBas,nBas,nspin),                                      &
             dipole_int_aa(nBas,nBas,ncart),dipole_int_bb(nBas,nBas,ncart),                  &
             ERI_aaaa(nBas,nBas,nBas,nBas),ERI_aabb(nBas,nBas,nBas,nBas),                    & 
             ERI_bbbb(nBas,nBas,nBas,nBas))
    if (doCAP) then
      allocate(complex_CAP_MO(nBas,nBas,nspin))
    else
      allocate(complex_CAP_MO(0,0,0))
    endif
    allocate(complex_ERI_aaaa(0,0,0,0))
    allocate(complex_ERI_aabb(0,0,0,0))
    allocate(complex_ERI_bbbb(0,0,0,0))
    allocate(complex_cHF(0,0,0))
    allocate(complex_dipole_int_aa(0,0,0))
    allocate(complex_dipole_int_bb(0,0,0)) 
  end if
  
  allocate(cHF(nBas,nBas,nspin),eHF(nBas,nspin))

  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)

  if(doMOM) then
    print *, "For MOM reference, only the following methods are available: GW(not implemented yet), RPA (not implemented yet)"
  end if

! For the FCIDUMP case, read two-body integrals

  if (readFCIDUMP) then 
    call read_fcidump_2body(nBas,ERI_AO)
  endif  

  call wall_time(end_int)
  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

!---------------------!
! Hartree-Fock module !
!---------------------!
  if(doUHF .and. .not. doMOM) then
    
    call wall_time(start_HF)
    call UHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EUHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UHF = ',t_HF,' seconds'
    write(*,*)

  end if

  if(doUHF .and. doMOM) then
    
    if(guess_type /= 6) then
      call wall_time(start_HF)
      call UHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
               nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EUHF,eHF,cHF,PHF,FHF)
      call wall_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UHF = ',t_HF,' seconds'
      write(*,*)
    end if

    call wall_time(start_HF)
    call MOM_UHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EUHF,eHF,cHF,PHF,FHF,mom_occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MOM-UHF = ',t_HF,' seconds'
    write(*,*)

  end if
  
  if(docUHF .and. .not. doMOM) then
   
    call wall_time(start_HF)
    call cUHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_EUHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for cUHF = ',t_HF,' seconds'
    write(*,*)
  
  end if

  if(docUHF .and. doMOM) then
   
    if(guess_type /= 6) then
      call wall_time(start_HF)
      call cUHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
               nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_EUHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF)
      call wall_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for cUHF = ',t_HF,' seconds'
      write(*,*)
    end if

    call wall_time(start_HF)
    call MOM_cUHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_EUHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF,mom_occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MOM-cUHF = ',t_HF,' seconds'
    write(*,*)

  end if
  


!------------------!
! Kohn-Sham module !
!------------------!

! if(doKS) then

!   call wall_time(start_KS)
!   write(*,*)
!   write(*,*) 'KS module has been disabled for now! Sorry.'
!   write(*,*)
!   call eDFT(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc,nBas,nC, & 
!             nO,nV,nR,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell,            &
!             max_ang_mom,min_exponent,max_exponent,S,T,V,Hc,X,ERI_AO,dipole_int_AO,EUHF,eHF,cHF,PHF,Vxc)
!   call wall_time(end_KS)

!   t_KS = end_KS - start_KS
!   write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for KS = ',t_KS,' seconds'
!   write(*,*)

! end if

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

  call wall_time(start_AOtoMO)

  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)

  ! Read and transform dipole-related integrals
  if(.not. docUHF) then
    do ixyz=1,ncart
      call AOtoMO(nBas,nBas,cHF(:,:,1),dipole_int_AO(:,:,ixyz),dipole_int_aa(:,:,ixyz))
      call AOtoMO(nBas,nBas,cHF(:,:,2),dipole_int_AO(:,:,ixyz),dipole_int_bb(:,:,ixyz))
    end do 
    
    ! 4-index transform for (aa|aa) block
    
    call AOtoMO_ERI_UHF(1,1,nBas,cHF,ERI_AO,ERI_aaaa)
    
    ! 4-index transform for (aa|bb) block
    
    call AOtoMO_ERI_UHF(1,2,nBas,cHF,ERI_AO,ERI_aabb)
    
    ! 4-index transform for (bb|bb) block
    
    call AOtoMO_ERI_UHF(2,2,nBas,cHF,ERI_AO,ERI_bbbb)
  
    call wall_time(end_AOtoMO)
  else
    do ixyz=1,ncart
      call complex_AOtoMO(nBas,nBas,complex_cHF(:,:,1),dipole_int_AO(:,:,ixyz),complex_dipole_int_aa(1,1,ixyz))
      call complex_AOtoMO(nBas,nBas,complex_cHF(:,:,2),dipole_int_AO(:,:,ixyz),complex_dipole_int_bb(1,1,ixyz))
    end do 
    
    ! 4-index transform for (aa|aa) block
    
    call complex_AOtoMO_ERI_UHF(1,1,nBas,complex_cHF,ERI_AO,complex_ERI_aaaa)
    
    ! 4-index transform for (aa|bb) block
    
    call complex_AOtoMO_ERI_UHF(1,2,nBas,complex_cHF,ERI_AO,complex_ERI_aabb)
    
    ! 4-index transform for (bb|bb) block
    
    call complex_AOtoMO_ERI_UHF(2,2,nBas,complex_cHF,ERI_AO,complex_ERI_bbbb)
    
    if(doCAP) then
      call complex_AOtoMO(nBas,nBas,complex_cHF(:,:,1),CAP_AO,complex_CAP_MO(:,:,1))
      call complex_AOtoMO(nBas,nBas,complex_cHF(:,:,2),CAP_AO,complex_CAP_MO(:,:,2))
    end if
  
  end if
  
  call wall_time(end_AOtoMO)
  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)

!-----------------------------------!
! Stability analysis of HF solution !
!-----------------------------------!

  nS(:) = (nO(:) - nC(:))*(nV(:) - nR(:))

  if(dostab) then

    call wall_time(start_stab)
    call UHF_stability(nBas,nC,nO,nV,nR,nS,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb)
    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

  if(dosearch) then

    call wall_time(start_stab)
    call UHF_search(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
                    nBas,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,      &
                    dipole_int_aa,dipole_int_bb,X,EUHF,eHF,cHF,PHF,FHF)

    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

!-----------------------!
! Moller-Plesset module !
!-----------------------!

  doMP = doMP2 .or. doMP3

  if(doMP) then

    call wall_time(start_MP)
    call UMP(dotest,doMP2,doMP3,reg_MP,doMOM,nCVS,mom_occupations,nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EUHF,eHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP = ',t_MP,' seconds'
    write(*,*)

  end if

!------------------------!
! Coupled-cluster module !
!------------------------!

  doCC = doCCD .or. dopCCD .or. doDCD .or. doCCSD .or. doCCSDT .or. &  
         dodrCCD .or. dorCCD .or. docrCCD .or. dolCCD

  if(doCC) then

    call wall_time(start_CC)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CC = ',t_CC,' seconds'
    write(*,*)

  end if

!----------------------------------!
! Configuration interaction module !
!----------------------------------!

  doCI = doCIS .or. doCID .or. doCISD .or. doFCI

  if(doCI) then

    call wall_time(start_CI)
    call UCI(dotest,doCIS,doCIS_D,doCID,doCISD,doFCI,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS, & 
             ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,EUHF,cHF,S)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CI = ',t_CI,' seconds'
    write(*,*)

  end if

!-----------------------------------!
! Random-phase approximation module !
!-----------------------------------!

  doRPA = dophRPA .or. dophRPAx .or. docrRPA .or. doppRPA
  CVS = any(nCVS>0) .or. any(FC>0) .or. doMOM
  
  if(doRPA .and. doUHF) then

    call wall_time(start_RPA)
    call URPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,CVS,  & 
              nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,&
              dipole_int_aa,dipole_int_bb,eHF,cHF,S,mom_occupations)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

  if(doRPA .and. docUHF) then

    call wall_time(start_RPA)
    call complex_URPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,CVS,  & 
              nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,complex_EUHF,complex_ERI_aaaa,complex_ERI_aabb,complex_ERI_bbbb,&
              complex_dipole_int_aa,complex_dipole_int_bb,complex_CAP_MO,complex_eHF,complex_cHF,S,mom_occupations)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if
!-------------------------!
! Green's function module !
!-------------------------!

  doGF = doG0F2 .or. doevGF2 .or. doqsGF2 .or. doG0F3 .or. doevGF3

  if(doGF) then

    call wall_time(start_GF)
    call UGF(dotest,doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,renorm_GF,maxSCF_GF,thresh_GF,max_diis_GF,    &
             dophBSE,doppBSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,lin_GF,eta_GF,reg_GF,               &
             nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb, &
             dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!-----------!
! GW module !
!-----------!

  doGW = doG0W0 .or. doevGW .or. doqsGW 

  if(doGW .and. .not. CVS) then
   ! After implementation of CVS GW variants remove condition on CVS 
    call wall_time(start_GW)
    call UGW(dotest,doG0W0,doevGW,doqsGW,maxSCF_GW,thresh_GW,max_diis_GW,                                         & 
             doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
             lin_GW,eta_GW,reg_GW,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,                        &  
             ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GW = ',t_GW,' seconds'
    write(*,*)

  end if
  
!-----------------!
! T-matrix module !
!-----------------!

  doGT = doG0T0pp .or. doevGTpp .or. doqsGTpp .or. doG0T0eh .or. doevGTeh .or. doqsGTeh

  if(doGT) then
    
    call wall_time(start_GT)
    call UGT(dotest,doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,maxSCF_GT,thresh_GT,max_diis_GT, &
             doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip,     &
             lin_GT,eta_GT,reg_GT,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,                            &
             ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------!
!     Parquet module     !
!------------------------!

  if(doevParquet) then
    call wall_time(start_Parquet)
!   call RParquet(max_it_macro,conv_one_body,max_it_micro,conv_two_body,   &
!        nOrb,nC,nO,nV,nR,nS, &
!        eHF,ERI_MO)            
    write(*,*) 'Unrestricted version of ev parquet not yet implemented. Sorry.'
    call wall_time(end_Parquet)
  
    t_Parquet = end_Parquet - start_Parquet
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Parquet module = ', t_Parquet, ' seconds'
    write(*,*)

 end if

 if(doqsParquet) then
    call wall_time(start_Parquet)
!   call RParquet(max_it_macro,conv_one_body,max_it_micro,conv_two_body,   &
!        nOrb,nC,nO,nV,nR,nS, &
!        eHF,ERI_MO)            
    write(*,*) 'Unrestricted version of qs parquet not yet implemented. Sorry.'
    call wall_time(end_Parquet)
  
    t_Parquet = end_Parquet - start_Parquet
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Parquet module = ', t_Parquet, ' seconds'
    write(*,*)

  end if

end subroutine
