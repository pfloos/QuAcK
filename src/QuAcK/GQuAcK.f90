subroutine GQuAcK(doGHF,dostab,dosearch,doMP2,doMP3,dophRPA,dophRPAx,doppRPA,                        &
                  doG0W0,doevGW,doqsGW,doG0F2,doevGF2,doqsGF2,                                       &
                  nNuc,nBas,nC,nO,nV,nR,ENuc,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,ERI_AO,              &
                  maxSCF_HF,max_diis_HF,thresh_HF,level_shift,guess_type,mix,reg_MP,                 &
                  TDA,maxSCF_GF,max_diis_GF,thresh_GF,lin_GF,reg_GF,eta_GF, &
                  maxSCF_GW,max_diis_GW,thresh_GW,TDA_W,lin_GW,reg_GW,eta_GW,                        &
                  dophBSE,dophBSE2,doppBSE,dBSE,dTDA,doACFDT,exchange_kernel,doXBS)

  implicit none
  include 'parameters.h'

  logical,intent(in)            :: doGHF
  logical,intent(in)            :: dostab
  logical,intent(in)            :: dosearch
  logical,intent(in)            :: doMP2
  logical,intent(in)            :: doMP3
  logical,intent(in)            :: dophRPA,dophRPAx,doppRPA
  logical,intent(in)            :: doG0F2,doevGF2,doqsGF2
  logical,intent(in)            :: doG0W0,doevGW,doqsGW

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
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

  logical,intent(in)            :: reg_MP

  logical,intent(in)            :: TDA

  integer,intent(in)            :: maxSCF_GF,max_diis_GF
  double precision,intent(in)   :: thresh_GF
  logical,intent(in)            :: lin_GF,reg_GF
  double precision,intent(in)   :: eta_GF

  integer,intent(in)            :: maxSCF_GW,max_diis_GW
  double precision,intent(in)   :: thresh_GW
  logical,intent(in)            :: TDA_W,lin_GW,reg_GW
  double precision,intent(in)   :: eta_GW

  logical,intent(in)            :: dophBSE,dophBSE2,doppBSE,dBSE,dTDA
  logical,intent(in)            :: doACFDT,exchange_kernel,doXBS

! Local variables

  logical                       :: doMP,doRPA,doGF,doGW
  
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: start_MP     ,end_MP       ,t_MP
  double precision              :: start_RPA    ,end_RPA      ,t_RPA
  double precision              :: start_GF     ,end_GF       ,t_GF
  double precision              :: start_GW     ,end_GW       ,t_GW

  double precision,allocatable  :: cHF(:,:),epsHF(:),PHF(:,:)
  double precision              :: EHF
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: ERI_tmp(:,:,:,:)
  double precision,allocatable  :: Ca(:,:),Cb(:,:)

  integer                       :: ixyz
  integer                       :: nBas2
  integer                       :: nS

  write(*,*)
  write(*,*) '*******************************'
  write(*,*) '* Generalized Branch of QuAcK *'
  write(*,*) '*******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nBas2 = 2*nBas

  allocate(cHF(nBas2,nBas2),epsHF(nBas2),PHF(nBas2,nBas2),   &
           dipole_int_MO(nBas2,nBas2,ncart),ERI_MO(nBas2,nBas2,nBas2,nBas2))

!---------------------!
! Hartree-Fock module !
!---------------------!

  if(doGHF) then

    call wall_time(start_HF)
    call GHF(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nBas2,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHF,epsHF,cHF,PHF)
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

  call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Ca,Ca,Ca,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Cb,Ca,Cb,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Ca,Cb,Ca,ERI_AO,ERI_tmp)
  ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

  call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Cb,Cb,Cb,ERI_AO,ERI_tmp)
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
    call GHF_stability(nBas2,nC,nO,nV,nR,nS,epsHF,ERI_MO)
    call wall_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

  if(dosearch) then

    call wall_time(start_stab)  
    call GHF_search(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                    nBas,nBas2,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHF,epsHF,cHF,PHF)
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
    call GMP(doMP2,doMP3,reg_MP,nBas2,nC,nO,nV,nR,ERI_MO,ENuc,EHF,epsHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP = ',t_MP,' seconds'
    write(*,*)

  end if

!-----------------------------------!
! Random-phase approximation module !
!-----------------------------------!

  doRPA = dophRPA .or. dophRPAx .or. doppRPA

  if(doRPA) then

    call wall_time(start_RPA)
    call GRPA(dophRPA,dophRPAx,doppRPA,TDA,doACFDT,exchange_kernel,nBas2,nC,nO,nV,nR,nS,ENuc,EHF, & 
              ERI_MO,dipole_int_MO,epsHF,cHF,S)
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
    call GGF(doG0F2,doevGF2,doqsGF2,maxSCF_GF,thresh_GF,max_diis_GF,dophBSE,doppBSE,TDA,dBSE,dTDA,lin_GF,eta_GF,reg_GF, &
             nNuc,ZNuc,rNuc,ENuc,nBas2,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,epsHF)
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
    call GGW(doG0W0,doevGW,doqsGW,maxSCF_GW,thresh_GW,max_diis_GW,doACFDT,exchange_kernel,doXBS,    & 
             dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,lin_GW,eta_GW,reg_GW,nNuc,ZNuc,rNuc,ENuc, & 
             nBas,nBas2,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,epsHF)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GW = ',t_GW,' seconds'
    write(*,*)

  end if

end subroutine
