subroutine GQuAcK(working_dir,dotest,doGHF,dostab,dosearch,doMP2,doMP3,doCCD,dopCCD,doDCD,doCCSD,doCCSDT, &
                  dodrCCD,dorCCD,docrCCD,dolCCD,dophRPA,dophRPAx,docrRPA,doppRPA,                         &
                  doG0W0,doevGW,doqsGW,doG0F2,doevGF2,doqsGF2,doG0T0pp,doevGTpp,doqsGTpp,doParquet,&
                  nNuc,nBas,nC,nO,nV,nR,ENuc,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,                          &
                  maxSCF_HF,max_diis_HF,thresh_HF,level_shift,guess_type,mix,reg_MP,                      &
                  maxSCF_CC,max_diis_CC,thresh_CC,                                                        &
                  TDA,maxSCF_GF,max_diis_GF,thresh_GF,lin_GF,reg_GF,eta_GF,                               &
                  maxSCF_GW,max_diis_GW,thresh_GW,TDA_W,lin_GW,reg_GW,eta_GW,                             &
                  maxSCF_GT,max_diis_GT,thresh_GT,TDA_T,lin_GT,reg_GT,eta_GT,                             &
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS,                       &
                  TDAeh,TDApp,max_diis_1b,max_diis_2b,max_it_1b,conv_1b,max_it_2b,conv_2b,lin_parquet,reg_parquet)

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doGHF
  logical,intent(in)            :: dostab
  logical,intent(in)            :: dosearch
  logical,intent(in)            :: doMP2
  logical,intent(in)            :: doMP3
  logical,intent(in)            :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical,intent(in)            :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical,intent(in)            :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical,intent(in)            :: doG0F2,doevGF2,doqsGF2
  logical,intent(in)            :: doG0W0,doevGW,doqsGW
  logical,intent(in)            :: doG0T0pp,doevGTpp,doqsGTpp
  logical,intent(in)            :: doParquet

  integer,intent(in)            :: nNuc,nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: ZNuc(nNuc),rNuc(nNuc,ncart)

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

  logical,intent(in)            :: reg_MP

  integer,intent(in)            :: maxSCF_CC,max_diis_CC
  double precision,intent(in)   :: thresh_CC

  logical,intent(in)            :: TDA

  integer,intent(in)            :: maxSCF_GF,max_diis_GF
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

  integer,intent(in)            :: max_it_1b,max_it_2b
  double precision,intent(in)   :: conv_1b,conv_2b
  integer,intent(in)            :: max_diis_1b,max_diis_2b
  logical,intent(in)            :: TDAeh,TDApp
  double precision,intent(in)   :: reg_parquet
  logical,intent(in)            :: lin_parquet

! Local variables

  logical                       :: doMP,doCC,doRPA,doGF,doGW,doGT
  
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: start_MP     ,end_MP       ,t_MP
  double precision              :: start_CC     ,end_CC       ,t_CC
  double precision              :: start_RPA    ,end_RPA      ,t_RPA
  double precision              :: start_GF     ,end_GF       ,t_GF
  double precision              :: start_GW     ,end_GW       ,t_GW
  double precision              :: start_GT     ,end_GT       ,t_GT
  double precision              :: start_Parquet,end_Parquet  ,t_Parquet

  double precision              :: start_int, end_int, t_int 
  double precision,allocatable  :: cHF(:,:),eHF(:),PHF(:,:),FHF(:,:)
  double precision              :: EGHF
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: ERI_tmp(:,:,:,:)
  double precision,allocatable  :: Ca(:,:),Cb(:,:)

  integer                       :: ixyz
  integer                       :: nBas2
  integer                       :: nS

  double precision,allocatable  :: eGW(:)

  write(*,*)
  write(*,*) '*******************************'
  write(*,*) '* Generalized Branch of QuAcK *'
  write(*,*) '*******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nBas2 = 2*nBas

  allocate(cHF(nBas2,nBas2),eHF(nBas2),PHF(nBas2,nBas2),FHF(nBas2,nBas2), &
           dipole_int_MO(nBas2,nBas2,ncart),ERI_MO(nBas2,nBas2,nBas2,nBas2))

  allocate(eGW(nBas2))

  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)
  call wall_time(end_int)
  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)


!---------------------!
! Hartree-Fock module !
!---------------------!

  if(doGHF) then

    call wall_time(start_HF)
    call GHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nBas2,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EGHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UHF = ',t_HF,' seconds'
    write(*,*)

  end if

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

  call wall_time(start_AOtoMO)

  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)

  allocate(Ca(nBas,nBas2),Cb(nBas,nBas2),ERI_tmp(nBas2,nBas2,nBas2,nBas2))

  Ca(:,:) = cHF(1:nBas,1:nBas2)
  Cb(:,:) = cHF(nBas+1:nBas2,1:nBas2)

  ! Transform dipole-related integrals

  do ixyz=1,ncart
    call AOtoMO_GHF(nBas,nBas2,Ca,Cb,dipole_int_AO(:,:,ixyz),dipole_int_MO(:,:,ixyz))
  end do 
  
  ! 4-index transform 

  call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Ca,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Cb,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Ca,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Cb,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

  deallocate(Ca,Cb,ERI_tmp)

  call wall_time(end_AOtoMO)

  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)

!-----------------------------------!
! Stability analysis of HF solution !
!-----------------------------------!

  nS = (nO - nC)*(nV - nR)

  if(dostab) then

    call wall_time(start_stab)
    call GHF_stability(nBas2,nC,nO,nV,nR,nS,eHF,ERI_MO)
    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

  if(dosearch) then

    call wall_time(start_stab)  
    call GHF_search(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                    nBas,nBas2,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,      & 
                    X,EGHF,eHF,cHF,PHF,FHF)
    call wall_time(end_stab)    

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

!-----------------------!
! Moller-Plesset module !
!-----------------------!

  doMP = doMP2 

  if(doMP) then

    call wall_time(start_MP)
    call GMP(dotest,doMP2,doMP3,reg_MP,nBas2,nC,nO,nV,nR,ERI_MO,ENuc,EGHF,eHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP = ',t_MP,' seconds'
    write(*,*)

  end if

!------------------------!
! Coupled-cluster module !
!------------------------!
    
  doCC = doCCD .or. doCCSD .or. doCCSDT .or. dodrCCD .or. dorCCD .or. docrCCD .or. dolCCD

  if(doCC) then

    call wall_time(start_CC)
    call GCC(dotest,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,dodrCCD,dorCCD,docrCCD,dolCCD, &
             maxSCF_CC,thresh_CC,max_diis_CC,nBas2,nC,nO,nV,nR,ERI_MO,ENuc,EGHF,eHF)
    call wall_time(end_CC)
  
    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CC = ',t_CC,' seconds'
    write(*,*)
    
  end if 

!-----------------------------------!
! Random-phase approximation module !
!-----------------------------------!

  doRPA = dophRPA .or. dophRPAx .or. docrRPA .or. doppRPA

  if(doRPA) then

    call wall_time(start_RPA)
    call GRPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!-------------------------!
! Green's function module !
!-------------------------!

  doGF = doG0F2 .or. doevGF2 .or. doqsGF2

  if(doGF) then

    call wall_time(start_GF)
    call GGF(dotest,doG0F2,doevGF2,doqsGF2,maxSCF_GF,thresh_GF,max_diis_GF,dophBSE,doppBSE,TDA,dBSE,dTDA,lin_GF,eta_GF,reg_GF, &
             nNuc,ZNuc,rNuc,ENuc,nBas2,nC,nO,nV,nR,nS,EGHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!-----------!
! GW module !
!-----------!

  doGW = doG0W0 .or. doevGW .or. doqsGW

  if(doGW) then
    
    call wall_time(start_GW)
    call GGW(dotest,doG0W0,doevGW,doqsGW,maxSCF_GW,thresh_GW,max_diis_GW,doACFDT,exchange_kernel,doXBS,       & 
             dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,lin_GW,eta_GW,reg_GW,nNuc,ZNuc,rNuc,ENuc,           & 
             nBas,nBas2,nC,nO,nV,nR,nS,EGHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF, &
             eGW)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GW = ',t_GW,' seconds'
    write(*,*)

 end if

!-----------------!
! T-matrix module !
!-----------------!

  doGT = doG0T0pp .or. doevGTpp .or. doqsGTpp

  if(doGT) then
    call wall_time(start_GT)
    call GGT(dotest,doG0T0pp,doevGTpp,doqsGTpp,  &
             maxSCF_GT,thresh_GT,max_diis_GT,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE, &
             TDA_T,TDA,dBSE,dTDA,lin_GT,eta_GT,reg_GT,nNuc,ZNuc,rNuc,ENuc,           &
             nBas,nBas2,nC,nO,nV,nR,nS,EGHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,                   &
             dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GT = ',t_GT,' seconds'
    write(*,*)

 end if

!------------------------!
!     Parquet module     !
!------------------------!

  if(doParquet) then
    call wall_time(start_Parquet)
    call GParquet(TDAeh,TDApp,max_diis_1b,max_diis_2b,lin_parquet,reg_parquet,ENuc,max_it_1b,conv_1b,max_it_2b,conv_2b, & 
                  nBas2,nC,nO,nV,nR,nS,EGHF,eHF,ERI_MO)
    call wall_time(end_Parquet)
  
    t_Parquet = end_Parquet - start_Parquet
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Parquet module = ',t_Parquet,' seconds'
    write(*,*)

  end if

end subroutine
