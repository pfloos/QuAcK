program QuAcK

  implicit none
  include 'parameters.h'

  logical                       :: unrestricted = .false.
  logical                       :: doHF,doRHF,doUHF,doRMOM,doUMOM
  logical                       :: dostab
  logical                       :: doKS
  logical                       :: doMP,doMP2,doMP3
  logical                       :: doCC,doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical                       :: dodrCCD,dorCCD,docrCCD,dolCCD
  logical                       :: doCI,doCIS,doCIS_D,doCID,doCISD,doFCI
  logical                       :: doRPA,dophRPA,dophRPAx,docrRPA,doppRPA
  logical                       :: doGF,doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3
  logical                       :: doGW,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,doSRGqsGW
  logical                       :: doGT,doG0T0pp,doevGTpp,doqsGTpp
  logical                       :: doG0T0eh,doevGTeh,doqsGTeh

  integer                       :: nNuc,nBas
  integer                       :: nEl(nspin)
  integer                       :: nC(nspin)
  integer                       :: nO(nspin)
  integer                       :: nV(nspin)
  integer                       :: nR(nspin)
  integer                       :: nS(nspin)
  double precision              :: ENuc,EHF

  double precision,allocatable  :: ZNuc(:),rNuc(:,:)
  double precision,allocatable  :: cHF(:,:,:),epsHF(:,:),PHF(:,:,:)

  logical                       :: doACFDT
  logical                       :: exchange_kernel
  logical                       :: doXBS

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: dipole_int_AO(:,:,:)
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: dipole_int_aa(:,:,:)
  double precision,allocatable  :: dipole_int_bb(:,:,:)
  double precision,allocatable  :: F_AO(:,:)
  double precision,allocatable  :: F_MO(:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  integer                       :: ixyz
  integer                       :: bra1,bra2
  integer                       :: ket1,ket2
  double precision,allocatable  :: ERI_MO_aaaa(:,:,:,:)
  double precision,allocatable  :: ERI_MO_aabb(:,:,:,:)
  double precision,allocatable  :: ERI_MO_bbbb(:,:,:,:)

  double precision              :: start_QuAcK  ,end_QuAcK    ,t_QuAcK
  double precision              :: start_int    ,end_int      ,t_int  
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_KS     ,end_KS       ,t_KS
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: start_MP     ,end_MP       ,t_MP 
  double precision              :: start_CC     ,end_CC       ,t_CC 
  double precision              :: start_CI     ,end_CI       ,t_CI 
  double precision              :: start_RPA    ,end_RPA      ,t_RPA 
  double precision              :: start_GF     ,end_GF       ,t_GF 
  double precision              :: start_GW     ,end_GW       ,t_GW 
  double precision              :: start_GT     ,end_GT       ,t_GT

  integer                       :: maxSCF_HF,max_diis_HF
  double precision              :: thresh_HF,level_shift
  logical                       :: DIIS_HF,guess_type,ortho_type,mix

  logical                       :: regMP

  integer                       :: maxSCF_CC,max_diis_CC
  double precision              :: thresh_CC
  logical                       :: DIIS_CC

  logical                       :: singlet
  logical                       :: triplet
  logical                       :: spin_conserved
  logical                       :: spin_flip
  logical                       :: TDA

  integer                       :: maxSCF_GF,max_diis_GF,renormGF
  double precision              :: thresh_GF
  logical                       :: DIIS_GF,linGF,regGF
  double precision              :: eta_GF

  integer                       :: maxSCF_GW,max_diis_GW
  double precision              :: thresh_GW
  logical                       :: DIIS_GW,TDA_W,linGW,regGW
  double precision              :: eta_GW

  integer                       :: maxSCF_GT,max_diis_GT
  double precision              :: thresh_GT
  logical                       :: DIIS_GT,TDA_T,linGT,regGT
  double precision              :: eta_GT

  logical                       :: dophBSE,dophBSE2,doppBSE,dBSE,dTDA


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

  call wall_time(start_QuAcK)

! Which calculations do you want to do?

  call read_methods(doRHF,doUHF,doRMOM,doUMOM,doKS,    &
                    doMP2,doMP3,                       &
                    doCCD,dopCCD,doDCD,doCCSD,doCCSDT, &
                    dodrCCD,dorCCD,docrCCD,dolCCD,     &
                    doCIS,doCIS_D,doCID,doCISD,doFCI,  & 
                    dophRPA,dophRPAx,docrRPA,doppRPA,  &
                    doG0F2,doevGF2,doqsGF2,            & 
                    doG0F3,doevGF3,                    &
                    doG0W0,doevGW,doqsGW,doSRGqsGW,    &
                    doufG0W0,doufGW,                   &
                    doG0T0pp,doevGTpp,doqsGTpp,        &
                    doG0T0eh,doevGTeh,doqsGTeh)

! Read options for methods

  call read_options(maxSCF_HF,thresh_HF,DIIS_HF,max_diis_HF,guess_type,ortho_type,mix,level_shift,dostab, &
                    regMP,                                                                              &
                    maxSCF_CC,thresh_CC,DIIS_CC,max_diis_CC,                                              &
                    TDA,singlet,triplet,spin_conserved,spin_flip,                                       &
                    maxSCF_GF,thresh_GF,DIIS_GF,max_diis_GF,linGF,eta_GF,renormGF,regGF,                  &
                    maxSCF_GW,thresh_GW,DIIS_GW,max_diis_GW,linGW,eta_GW,regGW,TDA_W,                     &  
                    maxSCF_GT,thresh_GT,DIIS_GT,max_diis_GT,linGT,eta_GT,regGT,TDA_T,                     & 
                    doACFDT,exchange_kernel,doXBS,                                                      &
                    dophBSE,dophBSE2,doppBSE,dBSE,dTDA)

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

  call read_molecule(nNuc,nEl,nO,nC,nR)
  allocate(ZNuc(nNuc),rNuc(nNuc,ncart))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

!------------------------------------------------------------------------
! Read basis set information from PySCF
!------------------------------------------------------------------------

  call read_basis_pyscf (nBas,nO,nV)
  nS(:) = (nO(:) - nC(:))*(nV(:) - nR(:))

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(cHF(nBas,nBas,nspin),epsHF(nBas,nspin),PHF(nBas,nBas,nspin),S(nBas,nBas),T(nBas,nBas),              &
           V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas),ERI_AO(nBas,nBas,nBas,nBas),dipole_int_AO(nBas,nBas,ncart), & 
           dipole_int_MO(nBas,nBas,ncart),F_AO(nBas,nBas))

! Read integrals

  call wall_time(start_int)

  call read_integrals(nBas,S,T,V,Hc,ERI_AO)
  call read_dipole_integrals(nBas,dipole_int_AO)

  call wall_time(end_int)

  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading integrals = ',t_int,' seconds'
  write(*,*)

! Compute orthogonalization matrix

  call orthogonalization_matrix(ortho_type,nBas,S,X)

!------------------------------------------------------------------------
! Hartree-Fock module
!------------------------------------------------------------------------

  doHF = doRHF .or. doUHF .or. doRMOM .or. doUMOM

  if(doHF) then

    call wall_time(start_HF)
    call HF(doRHF,doUHF,doRMOM,doUMOM,unrestricted,maxSCF_HF,thresh_HF,max_diis_HF, & 
            guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc,nBas,nO,S,T,V,Hc,F_AO, & 
            ERI_AO,dipole_int_AO,X,EHF,epsHF,cHF,PHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HF = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Kohn-Sham module
!------------------------------------------------------------------------

  if(doKS) then

    ! Switch on the unrestricted flag
    unrestricted = .true.

    call cpu_time(start_KS)
    write(*,*)
    write(*,*) 'KS module has been disabled for now! Sorry.'
    write(*,*)
!   call eDFT(maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc,nBas,nEl,nC, & 
!             nO,nV,nR,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell,            &
!             max_ang_mom,min_exponent,max_exponent,S,T,V,Hc,X,ERI_AO,dipole_int_AO,EHF,epsHF,cHF,PHF,Vxc)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for KS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! AO to MO integral transform for post-HF methods
!------------------------------------------------------------------------

  call wall_time(start_AOtoMO)

  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)

  if(unrestricted) then

    ! Read and transform dipole-related integrals
  
    allocate(dipole_int_aa(nBas,nBas,ncart),dipole_int_bb(nBas,nBas,ncart))
    dipole_int_aa(:,:,:) = dipole_int_AO(:,:,:)
    dipole_int_bb(:,:,:) = dipole_int_AO(:,:,:)
    do ixyz=1,ncart
        call AOtoMO_transform(nBas,cHF(:,:,1),dipole_int_aa(:,:,ixyz))
        call AOtoMO_transform(nBas,cHF(:,:,2),dipole_int_bb(:,:,ixyz))
    end do 

    ! Memory allocation
   
    allocate(ERI_MO_aaaa(nBas,nBas,nBas,nBas),ERI_MO_aabb(nBas,nBas,nBas,nBas),ERI_MO_bbbb(nBas,nBas,nBas,nBas))
   
    ! 4-index transform for (aa|aa) block
   
    bra1 = 1
    bra2 = 1
    ket1 = 1
    ket2 = 1
    call AOtoMO_integral_transform(bra1,bra2,ket1,ket2,nBas,cHF,ERI_AO,ERI_MO_aaaa)
    
    ! 4-index transform for (aa|bb) block
   
    bra1 = 1
    bra2 = 1
    ket1 = 2
    ket2 = 2
    call AOtoMO_integral_transform(bra1,bra2,ket1,ket2,nBas,cHF,ERI_AO,ERI_MO_aabb)
   
    ! 4-index transform for (bb|bb) block
   
    bra1 = 2
    bra2 = 2
    ket1 = 2
    ket2 = 2
    call AOtoMO_integral_transform(bra1,bra2,ket1,ket2,nBas,cHF,ERI_AO,ERI_MO_bbbb)

  else

    ! Memory allocation
   
    allocate(ERI_MO(nBas,nBas,nBas,nBas))
    allocate(F_MO(nBas,nBas))

    ! Read and transform dipole-related integrals
  
    dipole_int_MO(:,:,:) = dipole_int_AO(:,:,:)
    do ixyz=1,ncart
      call AOtoMO_transform(nBas,cHF,dipole_int_MO(:,:,ixyz))
    end do 

    ! 4-index transform 
   
    bra1 = 1
    bra2 = 1
    ket1 = 1
    ket2 = 1
    call AOtoMO_integral_transform(bra1,bra2,ket1,ket2,nBas,cHF,ERI_AO,ERI_MO)

    F_MO(:,:) = F_AO(:,:)
    call AOtoMO_transform(nBas,cHF,F_MO)

  end if

  call wall_time(end_AOtoMO)

  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)

!------------------------------------------------------------------------
! Stability analysis of HF solution
!------------------------------------------------------------------------

  if(dostab) then

    call cpu_time(start_stab)
    if(unrestricted) then
      call UHF_stability(nBas,nC,nO,nV,nR,nS,epsHF,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb)
    else
      call RHF_stability(nBas,nC,nO,nV,nR,nS,epsHF,ERI_MO)
    end if
    call cpu_time(end_stab)

    t_stab = end_stab - start_stab
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for stability analysis = ',t_stab,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Moller-Plesset module
!------------------------------------------------------------------------

  doMP = doMP2 .or. doMP3

  if(doMP) then

    call cpu_time(start_MP)
    call MP(doMP2,doMP3,unrestricted,regMP,nBas,nC,nO,nV,nR,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb,ENuc,EHF,epsHF)
    call cpu_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP = ',t_MP,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Coupled-cluster module
!------------------------------------------------------------------------

  doCC = doCCD .or. dopCCD .or. doDCD .or. doCCSD .or. doCCSDT .or. &  
         dodrCCD .or. dorCCD .or. docrCCD .or. dolCCD

  if(doCC) then

    call cpu_time(start_CC)
    call CC(doCCD,dopCCD,doDCD,doCCSD,doCCSDT,dodrCCD,dorCCD,docrCCD,dolCCD, & 
            maxSCF_CC,thresh_CC,max_diis_CC,nBas,nC,nO,nV,nR,ERI_MO,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CC = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Configuration interaction module
!------------------------------------------------------------------------

  doCI = doCIS .or. doCID .or. doCISD .or. doFCI

  if(doCI) then

    call cpu_time(start_CI)
    call CI(doCIS,doCIS_D,doCID,doCISD,doFCI,unrestricted,singlet,triplet,spin_conserved,spin_flip,                   &
            nBas,nC,nO,nV,nR,nS,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb,dipole_int_MO,dipole_int_aa,dipole_int_bb, &
            epsHF,EHF,cHF,S,F_MO)
    call cpu_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CI = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Random-phase approximation module
!------------------------------------------------------------------------

  doRPA = dophRPA .or. dophRPAx .or. docrRPA .or. doppRPA

  if(doRPA) then

    call cpu_time(start_RPA)
    call RPA(dophRPA,dophRPAx,docrRPA,doppRPA,unrestricted,                           & 
             TDA,doACFDT,exchange_kernel,singlet,triplet,spin_conserved,spin_flip,    & 
             nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb, & 
             dipole_int_MO,dipole_int_aa,dipole_int_bb,epsHF,cHF,S)
    call cpu_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Green's function module
!------------------------------------------------------------------------

  doGF = doG0F2 .or. doevGF2 .or. doqsGF2 .or. doG0F3 .or. doevGF3

  if(doGF) then

    call cpu_time(start_GF)
    call GF(doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,unrestricted,renormGF,maxSCF_GF,thresh_GF,max_diis_GF,              &
            dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,spin_conserved,spin_flip,linGF,eta_GF,regGF,                &
            nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb, &
            dipole_int_AO,dipole_int_MO,dipole_int_aa,dipole_int_bb,PHF,cHF,epsHF)
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! GW module
!------------------------------------------------------------------------

  doGW = doG0W0 .or. doevGW .or. doqsGW .or. doufG0W0 .or. doufGW .or. doSRGqsGW

  if(doGW) then
    
    call cpu_time(start_GW)
    call GW(doG0W0,doevGW,doqsGW,doufG0W0,doufGW,doSRGqsGW,unrestricted,maxSCF_GW,thresh_GW,max_diis_GW,doACFDT,         &
            exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,singlet,triplet,spin_conserved,spin_flip, &
            linGW,eta_GW,regGW,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,                                   &  
            ERI_AO,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb,dipole_int_AO,dipole_int_MO,dipole_int_aa,dipole_int_bb,   & 
            PHF,cHF,epsHF)
    call cpu_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GW = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! T-matrix module
!------------------------------------------------------------------------

  doGT = doG0T0pp .or. doevGTpp .or. doqsGTpp .or. doG0T0eh .or. doevGTeh .or. doqsGTeh

  if(doGT) then
    
    call cpu_time(start_GT)
    call GT(doG0T0pp,doevGTpp,doqsGTpp,doG0T0eh,doevGTeh,doqsGTeh,unrestricted,maxSCF_GT,thresh_GT,max_diis_GT,doACFDT,  &
            exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet,spin_conserved,spin_flip, &
            linGT,eta_GT,regGT,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,                                   &
            ERI_AO,ERI_MO,ERI_MO_aaaa,ERI_MO_aabb,ERI_MO_bbbb,dipole_int_AO,dipole_int_MO,dipole_int_aa,dipole_int_bb,   & 
            PHF,cHF,epsHF)
    call cpu_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! End of QuAcK
!------------------------------------------------------------------------

  call wall_time(end_QuAcK)

  t_QuAcK = end_QuAcK - start_QuAcK
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for QuAcK = ',t_QuAcK,' seconds'
  write(*,*)

end program 
