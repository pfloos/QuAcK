subroutine RQuAcK(working_dir,use_gpu,dotest,doRHF,doROHF,docRHF,doeRHF,doMOM,                                              &
                  dostab,dosearch,doaordm,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,                                    &
                  dodrCCD,dorCCD,docrCCD,dolCCD,doCIS,doCIS_D,doCID,doCISD,doFCI,dophRPA,dophRPAx,docrRPA,doppRPA,doOO,     & 
                  doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,doG0W0,doevGW,doqsGW,                                               &
                  doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp,doG0T0eh,doevGTeh,doqsGTeh,doevParquet,doqsParquet,                 & 
                  docG0W0,docG0F2,doscGW,doscGF2,                                                                           & 
                  doCAP,readFCIDUMP,restart_scGW,restart_scGF2,verbose_scGW,verbose_scGF2,chem_pot_scG,                     & 
                  do_IPEA_ADC2,do_IPEA_ADC3,do_SOSEX,do_2SOSEX,do_G3W2,                                                     & 
                  do_ADC_GW,do_ADC_2SOSEX,do_ADC3_G3W2,do_ADC3x_G3W2,do_ADC4_G3W2,                                          &
                  nNuc,nBas,nOrb,nC,nO,nV,nR,nCVS,FC,ENuc,ZNuc,rNuc,                                                        &
                  S,T,V,Hc,CAP_AO,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,eweight,eforward,             &
                  mom_occupations,writeMOs,                                                                                 &
                  guess_type,mix,reg_MP,maxSCF_CC,max_diis_CC,thresh_CC,singlet,triplet,TDA,                                &
                  maxIter_OO,thresh_OO,dRPA_OO,mu_OO,diagHess_OO,                                                           &
                  maxSCF_GF,max_diis_GF,renorm_GF,thresh_GF,lin_GF,reg_GF,eta_GF,maxSCF_GW,max_diis_GW,thresh_GW,           & 
                  TDA_W,lin_GW,reg_GW,eta_GW,do_linDM_GW,do_linDM_GF,                                                       &
                  maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,do_linDM_GT,                                   & 
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS,                                         &
                  TDAeh,TDApp,max_diis_1b,max_diis_2b,max_it_1b,conv_1b,max_it_2b,conv_2b,lin_parquet,reg_1b,reg_2b,reg_PA, &
                  nfreqs,ntimes,wcoord,wweight,                                                                             &
                  do_dyson,diag_approx,sig_inf,lin_ADC,reg_ADC,eta_ADC)      

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: use_gpu

  logical,intent(in)            :: dotest
  logical,intent(in)            :: readFCIDUMP

  logical,intent(in)            :: doRHF,doROHF,docRHF,doeRHF
  logical,intent(in)            :: doMOM,writeMOs
  logical,intent(in)            :: dostab
  logical,intent(in)            :: dosearch
  logical,intent(in)            :: doaordm
  logical,intent(in)            :: doscGW
  logical,intent(in)            :: doMP2,doMP3
  logical,intent(in)            :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical,intent(in)            :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical,intent(in)            :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical,intent(in)            :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical,intent(in)            :: doOO
  logical,intent(in)            :: doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,doscGF2
  logical,intent(inout)         :: doG0W0
  logical,intent(in)            :: doevGW,doqsGW
  logical,intent(in)            :: doG0T0pp,doevGTpp,doqsGTpp,doufG0T0pp
  logical,intent(in)            :: doG0T0eh,doevGTeh,doqsGTeh
  logical,intent(in)            :: docG0W0,docG0F2
  logical,intent(in)            :: doCAP
  logical,intent(in)            :: doevParquet,doqsParquet
  logical,intent(in)            :: do_IPEA_ADC2,do_IPEA_ADC3
  logical,intent(in)            :: do_SOSEX,do_2SOSEX,do_G3W2
  logical,intent(in)            :: do_ADC_GW,do_ADC_2SOSEX,do_ADC3_G3W2,do_ADC3x_G3W2,do_ADC4_G3W2

  integer,intent(in)            :: nNuc,nBas,nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  double precision,intent(inout):: ENuc
  integer,intent(in)            :: mom_occupations(nO,nspin)

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
  double precision,intent(in)   :: eweight
  logical,intent(in)            :: eforward

  logical,intent(in)            :: reg_MP

  integer,intent(in)            :: maxSCF_CC,max_diis_CC
  double precision,intent(in)   :: thresh_CC

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: TDA
 
  integer,intent(in)            :: maxIter_OO,thresh_OO,dRPA_OO,mu_OO,diagHess_OO

  integer,intent(in)            :: maxSCF_GF,max_diis_GF,renorm_GF
  double precision,intent(in)   :: thresh_GF
  logical,intent(in)            :: lin_GF,reg_GF
  double precision,intent(in)   :: eta_GF
  logical,intent(in)            :: do_linDM_GF
  logical,intent(in)            :: restart_scGF2
  logical,intent(in)            :: verbose_scGF2

  integer,intent(in)            :: maxSCF_GW,max_diis_GW
  double precision,intent(in)   :: thresh_GW
  logical,intent(in)            :: TDA_W,lin_GW,reg_GW
  double precision,intent(in)   :: eta_GW
  logical,intent(in)            :: do_linDM_GW
  logical,intent(in)            :: restart_scGW
  logical,intent(in)            :: verbose_scGW
  logical,intent(in)            :: chem_pot_scG

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
  
  logical,intent(in)            :: do_dyson,diag_approx,sig_inf,lin_ADC,reg_ADC
  double precision,intent(in)   :: eta_ADC

! Local variables

  logical                       :: doMP,doCC,doCI,doRPA,doGF,doGW,doGT,doADC
  logical                       :: file_exists
  logical                       :: no_fock
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
  double precision              :: start_ADC    ,end_ADC      ,t_ADC
  double precision              :: start_Parquet,end_Parquet  ,t_Parquet

  double precision              :: start_int, end_int, t_int
  double precision,allocatable  :: eHF(:)
  complex*16,allocatable        :: complex_eHF(:)
  double precision,allocatable  :: cHF(:,:)
  double precision,allocatable  :: cHF_tmp(:,:)
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
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: ixyz
  integer                       :: nS
  complex*16,allocatable        :: complex_ERI_MO(:,:,:,:)
  double precision,allocatable  :: eGW(:)

  logical                       :: doRDMs_numerically = .false.
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)

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

  allocate(cHF(nBas,nOrb))
  allocate(eGW(nOrb))
  
  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)

! For the FCIDUMP case, read two-body integrals

  if (readFCIDUMP) then 
    call read_fcidump_2body(nBas,ERI_AO)
  endif  
  
  call wall_time(end_int)
  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

  if(docRHF .and. (doRHF .or. doROHF)) then
    print *, "Complex Restricted-Hartree-Fock is not compatible with any other method than G0W0, evGW and qsGW !"
    stop
  end if
!---------------------!
! Hartree-Fock module !
!---------------------!

   if(doRHF .and. .not. doMOM) then

    call wall_time(start_HF)
    call RHF(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

  end if
  
  if(doRHF .and. doMOM) then

    if(guess_type /= 6) then
    
      call wall_time(start_HF)
      call RHF(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
               nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
      call wall_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
      write(*,*)

    end if

    call wall_time(start_HF)
    call MOM_RHF(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF,mom_occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

  end if

  if(doROHF .and. .not. doMOM) then

    call wall_time(start_HF)
    call ROHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
              nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ROHF = ',t_HF,' seconds'
    write(*,*)

  end if
  
  if(doROHF .and. doMOM) then
    if(guess_type /= 6) then
   
      call wall_time(start_HF)
      call ROHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
                nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
      call wall_time(end_HF)
    
    end if 
    
    call wall_time(start_HF)
    call MOM_ROHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
              nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF,mom_occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for MOM-ROHF = ',t_HF,' seconds'
    write(*,*)

  end if
  
  if(docRHF .and. .not. doMOM) then
  
    call wall_time(start_HF)
    call cRHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,ENuc, &
             nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_ERHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for cRHF = ',t_HF,' seconds'
    write(*,*)

  end if
  
  if(docRHF .and. doMOM) then
    if(guess_type /= 6) then 
      call wall_time(start_HF)
      call cRHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,ENuc, &
               nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_ERHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF)
      call wall_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for cRHF = ',t_HF,' seconds'
      write(*,*)
    endif
    call wall_time(start_HF)
    call MOM_cRHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,ENuc, &
             nBas,nO,S,T,V,ERI_AO,CAP_AO,X,complex_ERHF,complex_eHF,complex_cHF,complex_PHF,complex_FHF,mom_occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for cRHF = ',t_HF,' seconds'
    write(*,*)


  end if

!------------------------------!
! Ensemble Hartree-Fock module !
!------------------------------!

  if(doeRHF) then

    write(*,*)
    write(*,'(A)') ' Warning! The HF density and orbitals will be replaced by the ensemble ones.'
    write(*,*)
    ! Do an ensemble eRHF calculation 
    call wall_time(start_HF)
    call ensembleRHF(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                     nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eweight,eforward,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for eRHF = ',t_HF,' seconds'
    write(*,*)
 
  end if

!--------------!
! scGF2 module !
!--------------!

  if(doscGF2 .and. .not.docRHF) then

   allocate(vMAT(nBas*nBas,nBas*nBas))
   allocate(cHF_tmp(nBas,nOrb))
   cHF_tmp=cHF
   do iorb=1,nBas
    do jorb=1,nBas
     do korb=1,nBas
      do lorb=1,nBas
       vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_AO(iorb,jorb,korb,lorb)
      enddo
     enddo
    enddo
   enddo
   no_fock=.false.
   call scGF2_AO_itau_iw(nBas,nOrb,nO,maxSCF_GF,max_diis_GF,do_linDM_GF,restart_scGF2,verbose_scGF2,chem_pot_scG,no_fock, &
                         ENuc,Hc,S,PHF,cHF_tmp,eHF,nfreqs,wcoord,wweight,vMAT,ERI_AO)
   deallocate(vMAT,cHF_tmp)

  endif

!-------------!
! scGW module !
!-------------!

  if(doscGW .and. .not.docRHF) then

   allocate(vMAT(nBas*nBas,nBas*nBas))
   allocate(cHF_tmp(nBas,nOrb))
   cHF_tmp=cHF
   do iorb=1,nBas
    do jorb=1,nBas
     do korb=1,nBas
      do lorb=1,nBas
       vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_AO(iorb,jorb,korb,lorb)
      enddo
     enddo
    enddo
   enddo
   no_fock=.false.
!  Use scGHF to check the procedure and the convergence when Sigma_Hxc = Sigma_Hx (NOTE: Modify RHF to enforce Hcore or RH)
!   call scGHF_AO_itau_iw(nBas,nOrb,nO,maxSCF_GW,max_diis_GW,verbose_scGW,restart_scGW,chem_pot_scG, &
!                         ENuc,Hc,S,X,PHF,cHF_tmp,eHF,nfreqs,wcoord,wweight,vMAT)
   call scGW_AO_itau_iw(nBas,nOrb,nO,maxSCF_GW,max_diis_GW,do_linDM_GW,restart_scGW,verbose_scGW,chem_pot_scG,no_fock, &
                        ENuc,Hc,S,X,PHF,cHF_tmp,eHF,nfreqs,wcoord,wweight,vMAT,ERI_AO)
   deallocate(vMAT,cHF_tmp)

  endif


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
    call RHF_search(maxSCF_HF,doaordm,thresh_HF,max_diis_HF,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
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
  CVS = any(nCVS>0) .or. any(FC > 0) .or. doMOM

  if(doRPA .and. doRHF) then

    call wall_time(start_RPA)
    call RRPA(use_gpu,dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet,CVS, &
              nOrb,nC,nO,nV,nR,nS,nCVS(1),FC(1),ENuc,ERHF,ERI_MO,dipole_int_MO,eHF,mom_occupations(:,1))
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if
  
  if(doRPA .and. docRHF) then

    call wall_time(start_RPA)
    call complex_RRPA(use_gpu,dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet,CVS, &
              nOrb,nC,nO,nV,nR,nS,nCVS(1),FC(1),ENuc,complex_ERHF,complex_ERI_MO,complex_dipole_int_MO,complex_CAP_MO,complex_eHF,mom_occupations(:,1))
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if
!-------------------------------------------!
! Orbital optimisation module RPA/GW module !
!-------------------------------------------!
  if(doG0W0 .and. doOO) then

    call wall_time(start_GW)
    call OORG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
         lin_GW,eta_GW,reg_GW,nBas,nOrb,nC,nO,nV,nR,nS,                                                                &
         maxIter_OO,thresh_OO,dRPA_OO,mu_OO,diagHess_OO,                                                            &
         ENuc,ERHF,ERI_AO,ERI_MO,dipole_int_MO,eHF,cHF,                                                             &
         S,X,T,V,Hc,PHF,FHF,eGW) 
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Orbital optimisation and G0W0 = ',t_GW,' seconds'
    write(*,*)
    doG0W0 = .false.  
  
  end if
!-------------------------!
! Green's function module !
!-------------------------!

doGF = doG0F2 .or. doevGF2 .or. doqsGF2 .or. doG0F3 .or. doevGF3 .or. docG0F2

  if(doGF .and. .not. docRHF) then
    call wall_time(start_GF)
    call RGF(dotest,doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,renorm_GF,maxSCF_GF,           &
             thresh_GF,max_diis_GF,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,lin_GF, &
             eta_GF,reg_GF,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,complex_ERHF,    &
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

  doGW = doG0W0 .or. doevGW .or. doqsGW .or. docG0W0 

  if(doGW .and. .not. docRHF) then
    
    call wall_time(start_GW)
    call RGW(dotest,doG0W0,doevGW,doqsGW,maxSCF_GW,thresh_GW,max_diis_GW,                                & 
             doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,singlet,triplet, &
             lin_GW,eta_GW,reg_GW,do_linDM_GW,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,   &
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
             TDA_T,TDA,dBSE,dTDA,singlet,triplet,lin_GT,eta_GT,reg_GT,do_linDM_GT,nNuc,ZNuc,rNuc,ENuc,&
             nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,                    &
             dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------!
! ADC module !
!------------!

  doADC = do_IPEA_ADC2 .or. do_IPEA_ADC3 .or.       &
          do_SOSEX .or. do_2SOSEX .or. do_G3W2 .or. &
          do_ADC_GW .or. do_ADC_2SOSEX .or. do_ADC3_G3W2 .or. do_ADC3x_G3W2 .or. do_ADC4_G3W2

  if(doADC) then

    call wall_time(start_ADC)
    call R_ADC(dotest,                                               &
               do_IPEA_ADC2,do_IPEA_ADC3,                            & 
               do_SOSEX,do_2SOSEX,do_G3W2,                           & 
               do_ADC_GW,do_ADC_2SOSEX,                              &
               do_ADC3_G3W2,do_ADC3x_G3W2,do_ADC4_G3W2,              &
               TDA_W,TDA,singlet,triplet,lin_ADC,eta_ADC,reg_ADC,    &
               do_dyson,diag_approx,sig_inf,                         &
               nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,         &
               S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO, &
               ERHF,PHF,FHF,cHF,eHF)
    call wall_time(end_ADC)

    t_ADC = end_ADC - start_ADC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC = ',t_ADC,' seconds'
    write(*,*)

  end if

!----------------!
! Parquet module !
!----------------!

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

  if(doRDMs_numerically) then
    call R_rdm1_numerical(dotest,doaordm,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,eHF,cHF,PHF,FHF)    
  endif

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
