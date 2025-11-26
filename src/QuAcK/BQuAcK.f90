subroutine BQuAcK(working_dir,dotest,doaordm,doRHFB,doBRPA,dophRPA,doG0W0,doqsGW,doscGW,readFCIDUMP,nNuc,nBas,nOrb, &
                  nO,ENuc,eta,shift,restart_scGW,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,maxSCF,max_diis,thresh,         &
                  level_shift,guess_type,maxSCF_GW,max_diis_GW,thresh_GW,dolinGW,temperature,sigma,                 &
                  chem_pot_hf,restart_hfb,nfreqs,ntimes,wcoord,wweight,error_P,verbose_scGW,chem_pot_scG)

! Restricted branch of Bogoliubov QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in)  :: working_dir
                                 
  logical,intent(in)             :: dotest
  logical,intent(in)             :: doaordm
  logical,intent(in)             :: readFCIDUMP
  logical,intent(in)             :: error_P
  logical,intent(in)             :: verbose_scGW
  logical,intent(in)             :: chem_pot_scG
                                 
  logical,intent(in)             :: doRHFB
  logical,intent(in)             :: doBRPA
  logical,intent(in)             :: dophRPA
  logical,intent(in)             :: doG0W0
  logical,intent(in)             :: doqsGW
  logical,intent(in)             :: dolinGW
  logical,intent(in)             :: doscGW
  logical,intent(in)             :: restart_scGW

  logical,intent(in)             :: restart_hfb
  logical,intent(in)             :: chem_pot_hf
  integer,intent(in)             :: nNuc,nBas,nOrb
  integer,intent(in)             :: nO
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: ntimes
  double precision,intent(inout) :: ENuc
  double precision,intent(in)    :: eta
  double precision,intent(in)    :: shift
  double precision,intent(in)    :: temperature,sigma

  double precision,intent(in)    :: ZNuc(nNuc),rNuc(nNuc,ncart)
  double precision,intent(in)    :: wcoord(nfreqs),wweight(nfreqs)

  double precision,intent(inout) :: S(nBas,nBas)
  double precision,intent(inout) :: T(nBas,nBas)
  double precision,intent(inout) :: V(nBas,nBas)
  double precision,intent(inout) :: Hc(nBas,nBas)
  double precision,intent(inout) :: X(nBas,nOrb)
  double precision,intent(inout) :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)             :: maxSCF,max_diis
  integer,intent(in)             :: maxSCF_GW,max_diis_GW
  integer,intent(in)             :: guess_type
  double precision,intent(in)    :: thresh,level_shift
  double precision,intent(in)    :: thresh_GW

! Local variables

  logical                        :: no_fock
  logical                        :: file_exists

  integer                        :: verbose
  integer                        :: nOrb_twice
  integer                        :: ixyz
  integer                        :: iorb,jorb,korb,lorb
                                
  double precision               :: chem_pot,Val
  double precision               :: start_HF     ,end_HF       ,t_HF
  double precision               :: start_Ecorr  ,end_Ecorr    ,t_Ecorr
  double precision               :: start_qsGWB  ,end_qsGWB    ,t_qsGWB
  double precision               :: start_int    ,end_int      ,t_int
                                
  double precision,allocatable   :: eHF(:)
  double precision,allocatable   :: eQP_state(:)
  double precision,allocatable   :: U_QP(:,:)
  double precision,allocatable   :: MOCoef(:,:)
  double precision,allocatable   :: pMAT(:,:)
  double precision,allocatable   :: pMATcorr(:,:)
  double precision,allocatable   :: panomMAT(:,:)
  double precision,allocatable   :: panomMATcorr(:,:)
  double precision,allocatable   :: Fock(:,:)
  double precision,allocatable   :: Delta(:,:)
  double precision,allocatable   :: vMAT(:,:)
  double precision               :: EeleSD,Eelec,EcRPA,EcGM
  double precision,allocatable   :: dipole_int_MO(:,:,:)
  double precision,allocatable   :: ERI_AO(:,:,:,:)
  double precision,allocatable   :: ERI_MO(:,:,:,:)
                                
!

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  verbose=0
  nOrb_twice=nOrb+nOrb

  allocate(eHF(nOrb))

  allocate(MOCoef(nBas,nOrb))

  allocate(pMAT(nBas,nBas))
  allocate(panomMAT(nBas,nBas))
  allocate(Fock(nBas,nBas))
  allocate(Delta(nBas,nBas))

  allocate(eQP_state(nOrb_twice))
  allocate(U_QP(nOrb_twice,nOrb_twice))

  allocate(ERI_AO(nBas,nBas,nBas,nBas))

  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)

! For the FCIDUMP read two-body parameters integrals

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
   close(314)
  endif
  
  call wall_time(end_int)
  t_int = end_int - start_int

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

!--------------------------------!
! Hartree-Fock Bogoliubov module !
!--------------------------------!

  if(doRHFB) then

    ! Run first a RHF calculation 
    call wall_time(start_HF)
    call RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)
    
    ! Compute G0W0 AND/OR EcRPA and EcGM energies for RHF and qsGW
    if(dophRPA .or. doG0W0) then

     call wall_time(start_Ecorr)
     allocate(vMAT(nOrb*nOrb,nOrb*nOrb))
     allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
     call AOtoMO_ERI_RHF(nBas,nOrb,MOCoef,ERI_AO,ERI_MO)
     do iorb=1,nOrb
      do jorb=1,nOrb
       do korb=1,nOrb
        do lorb=1,nOrb
         vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
        enddo
       enddo
      enddo
     enddo
     deallocate(ERI_MO)
     if(dophRPA) then
      call EcRPA_EcGM_w_RHF(nOrb,nO,1,eHF,nfreqs,ntimes,wweight,wcoord,vMAT,EeleSD+ENuc, &
                            EcRPA,EcGM)
     endif
     ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr
     if(dolinGW .and. dophRPA) then
      allocate(pMATcorr(nBas,nBas))
      call linDyson_GW_RHF(nBas,nOrb,nO,MOCoef,eHF,nfreqs,wweight,wcoord,ERI_AO,vMAT, &
                          Enuc,EcGM,T,V,S,pMAT,pMATcorr)
      deallocate(pMATcorr)
     endif
     deallocate(vMAT)
     call wall_time(end_Ecorr)

     t_Ecorr = end_Ecorr - start_Ecorr
     write(*,*)
     write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Ecorr = ',t_Ecorr,' seconds'
     write(*,*)

    endif

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call RHFB(dotest,doaordm,maxSCF,thresh,max_diis,error_P,nNuc,ZNuc,rNuc,ENuc,          &
             nBas,nOrb,nOrb_twice,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,Eelec,eHF,MOCoef,    & 
             pMAT,panomMAT,Fock,Delta,temperature,sigma,chem_pot_hf,chem_pot,restart_hfb, &
             U_QP,eQP_state)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

  ! Compute G0W0 AND/OR EcRPA and EcGM energies for qsGWB or HFB
  if(dophRPA .or. doG0W0) then

   call wall_time(start_Ecorr)
   allocate(vMAT(nOrb*nOrb,nOrb*nOrb))
   allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
   call AOtoMO_ERI_RHF(nBas,nOrb,MOCoef,ERI_AO,ERI_MO)
   do iorb=1,nOrb
    do jorb=1,nOrb
     do korb=1,nOrb
      do lorb=1,nOrb
       vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
      enddo
     enddo
    enddo
   enddo
   deallocate(ERI_MO)
   if(dophRPA) then
    call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                           U_QP,Eelec+ENuc,EcRPA,EcGM)
   endif
   ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr and Panomcorr
   if(dolinGW .and. dophRPA) then
    allocate(pMATcorr(nBas,nBas),panomMATcorr(nBas,nBas))
    call linDyson_GW_RHFB(nBas,nOrb,nOrb_twice,MOCoef,eQP_state,nfreqs,wweight,wcoord,ERI_AO,vMAT,U_QP,&
                          Enuc,EcGM,sigma,T,V,S,X,pMAT,panomMAT,pMATcorr,panomMATcorr)
    deallocate(pMATcorr,panomMATcorr)
   endif
   deallocate(vMAT)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
   write(*,*)
   write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Ecorr = ',t_Ecorr,' seconds'
   write(*,*)

  endif

  ! scGW Bogoliubov
  if(doscGW) then
   allocate(vMAT(nBas*nBas,nBas*nBas))
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
   call scGWB_AO_itau_iw(nBas,nOrb,nOrb_twice,maxSCF_GW,max_diis_GW,dolinGW,restart_scGW,verbose_scGW,chem_pot_scG,no_fock, &
                         ENuc,Hc,S,X,pMAT,panomMAT,MOCoef,eQP_state,chem_pot,sigma,nfreqs,wcoord,wweight,U_QP,vMAT,ERI_AO)
   deallocate(vMAT)
  endif

!-----------------------------------------------------------------!
! Compute the EcRPA Bogoliubov using the FCIDUMP and occ. numbers !
!-----------------------------------------------------------------!

  if(doBRPA) then

   call wall_time(start_Ecorr)
   deallocate(pMAT,panomMAT,Fock,Delta)
   allocate(pMAT(nOrb,nOrb),panomMAT(nOrb,nOrb))
   allocate(Fock(nOrb,nOrb),Delta(nOrb,nOrb))
   allocate(vMAT(nOrb*nOrb,nOrb*nOrb),ERI_MO(nOrb,nOrb,nOrb,nOrb))
   call EcRPA_Bogoliubov_FCIDUMP(nO,nOrb,nOrb_twice,sigma,ERI_MO,vMAT,Fock,Delta,pMAT,panomMAT,eQP_state,U_QP, &
                                 chem_pot,ntimes,nfreqs,wcoord,wweight)
   deallocate(vMAT,ERI_MO)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
   write(*,*)
   write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Ecorr = ',t_Ecorr,' seconds'
   write(*,*)

  endif

! Memory deallocation
    
  deallocate(eHF)
  deallocate(MOCoef)
  deallocate(pMAT)
  deallocate(panomMAT)
  deallocate(Fock)
  deallocate(Delta)
  deallocate(eQP_state)
  deallocate(U_QP)
  deallocate(ERI_AO)

end subroutine
