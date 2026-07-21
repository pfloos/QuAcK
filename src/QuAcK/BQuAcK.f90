subroutine BQuAcK(working_dir,dotest,doaordm,doRHFB,doBRPA,dophRPA,dophRPAx,doMP2,doscGF2,doscGW,readFCIDUMP,nNuc,nBas,  &
                  nOrb,nO,ENuc,eta,shift,restart_scGW,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,maxSCF,max_diis,doscGHF,thresh, &
                  level_shift,guess_type,mix,TDA,maxSCF_GW,max_diis_GW,thresh_GW,dolinGW,dosign_XoB,temperature,sigma,   &
                  maxSCF_GF,max_diis_GF,thresh_GF,restart_scGF2,verbose_scGF2,                                           &
                  chem_pot_hf,restart_hfb,nfreqs,ntimes,wcoord,wweight,error_P,verbose_scGW,chem_pot_scG,writeMOs)

! Restricted branch of Bogoliubov QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in)  :: working_dir
                                 
  logical,intent(in)             :: dotest
  logical,intent(in)             :: doaordm
  logical,intent(in)             :: readFCIDUMP
  logical,intent(in)             :: error_P
  logical,intent(in)             :: verbose_scGW
  logical,intent(in)             :: verbose_scGF2
  logical,intent(in)             :: chem_pot_scG
  logical,intent(in)             :: TDA
                                 
  logical,intent(in)             :: doRHFB
  logical,intent(in)             :: doBRPA
  logical,intent(in)             :: dophRPAx
  logical,intent(in)             :: dophRPA
  logical,intent(in)             :: doMP2
  logical,intent(in)             :: dolinGW
  logical,intent(in)             :: dosign_XoB
  logical,intent(in)             :: doscGF2
  logical,intent(in)             :: doscGW
  logical,intent(in)             :: doscGHF
  logical,intent(in)             :: restart_scGW
  logical,intent(in)             :: restart_scGF2
  logical,intent(in)             :: writeMOs

  logical,intent(in)             :: restart_hfb
  logical,intent(in)             :: chem_pot_hf
  integer,intent(in)             :: nNuc,nBas,nOrb
  integer,intent(in)             :: nO
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: ntimes
  double precision,intent(inout) :: ENuc
  double precision,intent(in)    :: eta
  double precision,intent(in)    :: mix
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
  integer,intent(in)             :: maxSCF_GF,max_diis_GF
  double precision,intent(in)    :: thresh_GF

! Local variables

  logical                        :: no_fock
  logical                        :: file_exists

  integer                        :: verbose
  integer                        :: nOrb_twice
  integer                        :: ixyz
  integer                        :: iorb,jorb,korb,lorb
                                
  double precision               :: sign_XoB
  double precision               :: EcRPAx
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
  double precision               :: EeleSD,Eelec,EcRPA,EcGM,EcMP2
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

  sign_XoB=-1d0
  if(dosign_XoB) sign_XoB=1d0
  verbose=0
  nOrb_twice=nOrb+nOrb
  Eelec=0d0

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

  if (readFCIDUMP) then 
    call read_fcidump_2body(nBas,ERI_AO)
  endif  
  
  call wall_time(end_int)
  t_int = end_int - start_int

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

!--------------------------------!
! Hartree-Fock Bogoliubov module !
!--------------------------------!

  ! Test scGHF to do RHF after restricted Hartree
  if(doscGHF) then
   call RH(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
            nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
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
   call scGHF_AO_itau_iw(nBas,nOrb,nO,maxSCF,max_diis,0,restart_scGW,chem_pot_scG, &
                         ENuc,Hc,S,X,pMAT,MOCoef,eHF,nfreqs,wcoord,wweight,vMAT)
   deallocate(vMAT)
  endif

  if(doRHFB) then

    ! Run first a RHF calculation 
    call wall_time(start_HF)
    call RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

    ! Compute EcRPA and EcGM energies and lin-G for RHF
    if(dophRPA) then

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
     call EcRPA_EcGM_w_RHF(nOrb,nO,1,eHF,nfreqs,ntimes,wweight,wcoord,vMAT,EeleSD+ENuc, &
                           EcRPA,EcGM)
     ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr
     if(dolinGW) then
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

  ! Compute EcMP2 for RHFB
  if(doMP2) then
   call wall_time(start_Ecorr)
   call BMP2(nBas,nOrb,MOCoef,Hc,S,ERI_AO,chem_pot,sigma,U_QP,Eelec+ENuc,EcMP2)
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
   ! Direct MP2 term from Xo(i w)
   if(nfreqs>1) then
    call EcMP2_w_RHFB(nOrb,nOrb_twice,1,sign_XoB,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT,&
                      U_QP,Eelec+ENuc,EcMP2)
   endif
   deallocate(vMAT)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
   write(*,*)
   write(*,'(A65,1X,F9.3,A8)') 'Total wall time for EBMP2 = ',t_Ecorr,' seconds'
   write(*,*)
  endif

  ! scGF2 Bogoliubov
  if(doscGF2) then
   no_fock=.false.
   ! MRM: TODO we should fix the restricted version of Bog. scGF2 because we are using the generalized one
   call faked_Gen_scGF2B(nBas,nOrb,nOrb_twice,nfreqs,no_fock,Enuc,sigma,chem_pot,maxSCF_GF,max_diis_GF,thresh_GF,restart_scGF2,verbose_scGF2,chem_pot_scG, &
                         wcoord,wweight,Fock,Delta,Hc,S,X,pMAT,panomMAT,MOCoef,ERI_AO)
! TODO fix the spin adaptation
!   call scGF2B_AO_itau_iw(nBas,nOrb,nOrb_twice,maxSCF_GF,thresh_GF,max_diis_GF,restart_scGF2,verbose_scGF2,       &
!                          chem_pot_scG,no_fock,ENuc,Hc,S,X,pMAT,panomMAT,MOCoef,eQP_state,chem_pot,sigma,nfreqs,  &
!                          wcoord,wweight,U_QP,ERI_AO)
  endif

  ! Compute EcRPAx for RHFB
  if(dophRPAx .and. .false.) then
   call wall_time(start_Ecorr)
   call BRPAx(nBas,nOrb,TDA,MOCoef,Hc,S,ERI_AO,chem_pot,sigma,U_QP,Eelec+ENuc,EcRPAx)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
   write(*,*)
   write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ERPAx = ',t_Ecorr,' seconds'
   write(*,*)
  endif

  ! Compute EcRPA and EcGM energies and lin-G for RHFB
  if(dophRPA) then

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
   call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,sign_XoB,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                          U_QP,Eelec+ENuc,EcRPA,EcGM)
   ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr and Panomcorr
   if(dolinGW) then
    allocate(pMATcorr(nBas,nBas),panomMATcorr(nBas,nBas))
    call linDyson_GW_RHFB(nBas,nOrb,nOrb_twice,sign_XoB,MOCoef,eQP_state,nfreqs,wweight,wcoord,ERI_AO,vMAT,U_QP,&
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
   call scGWB_AO_itau_iw(nBas,nOrb,nOrb_twice,maxSCF_GW,thresh_GW,max_diis_GW,dolinGW,dophRPA,restart_scGW,verbose_scGW, &
                         chem_pot_scG,no_fock,ENuc,Hc,S,X,pMAT,panomMAT,MOCoef,eQP_state,chem_pot,sigma,sign_XoB,nfreqs, &
                         wcoord,wweight,U_QP,vMAT,ERI_AO)
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
   call EcRPA_Bogoliubov_FCIDUMP(nO,nOrb,nOrb_twice,sigma,sign_XoB,ERI_MO,vMAT,Fock,Delta,pMAT,panomMAT,eQP_state,U_QP, &
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

subroutine faked_Gen_scGF2B(nBas,nOrb,nOrb_twice,nfreqs,no_fock,Enuc,sigma,chem_pot,maxSCF_GF,max_diis_GF,thresh_GF,restart_scGF2,verbose_scGF2,chem_pot_scG, &
                            wcoord,wweight,Fock,Delta,Hc,S,X,pMAT,panomMAT,MOCoef,ERI_AO)

  logical,intent(in)             :: no_fock
  logical,intent(in)             :: restart_scGF2
  logical,intent(in)             :: verbose_scGF2
  logical,intent(in)             :: chem_pot_scG
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: maxSCF_GF,max_diis_GF
  integer,intent(in)             :: nBas,nOrb,nOrb_twice
  double precision,intent(in)    :: ENuc
  double precision,intent(in)    :: sigma
  double precision,intent(in)    :: chem_pot
  double precision,intent(in)    :: thresh_GF
  double precision,intent(in)    :: wcoord(nfreqs)
  double precision,intent(in)    :: wweight(nfreqs)
  double precision,intent(in)    :: S(nBas,nBas)
  double precision,intent(in)    :: X(nBas,nOrb)
  double precision,intent(in)    :: Hc(nBas,nBas)
  double precision,intent(inout) :: Fock(nBas,nBas)
  double precision,intent(inout) :: Delta(nBas,nBas)
  double precision,intent(in)    :: pMAT(nBas,nBas)
  double precision,intent(in)    :: panomMAT(nBas,nBas)
  double precision,intent(in)    :: MOCoef(nBas,nOrb)
  double precision,intent(in)    :: ERI_AO(nBas,nBas,nBas,nBas)

  integer                        :: ibas,jbas,kbas,lbas
  integer                        :: nBas_twice
  integer                        :: nBas_3times
  integer                        :: nBas_4times
  integer                        :: nOrb_3times
  integer                        :: nOrb_4times
  double precision               :: trace_1rdm
  double precision,allocatable   :: Gen_eQP_state(:)
  double precision,allocatable   :: J(:,:)
  double precision,allocatable   :: K(:,:)
  double precision,allocatable   :: Gen_Hc(:,:)
  double precision,allocatable   :: Gen_S(:,:)
  double precision,allocatable   :: Gen_R(:,:)
  double precision,allocatable   :: Gen_R_mo(:,:)
  double precision,allocatable   :: Gen_MOCoef(:,:)
  double precision,allocatable   :: Gen_H_HFB_ao(:,:)
  double precision,allocatable   :: Gen_MOCoef4(:,:)
  double precision,allocatable   :: Gen_U_QP(:,:)
  double precision,allocatable   :: sw_ERI_AO(:,:,:,:)
  double precision,allocatable   :: db_ERI_AO(:,:,:,:)

  nBas_twice=2*nBas
  nBas_3times=nBas+nBas_twice
  nBas_4times=2*nBas_twice
  nOrb_3times=nOrb+nOrb_twice
  nOrb_4times=2*nOrb_twice
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Gen_eQP_state(nOrb_4times))
  allocate(Gen_Hc(nBas_twice,nBas_twice))
  allocate(Gen_S(nBas_twice,nBas_twice))
  allocate(Gen_R(nBas_4times,nBas_4times))
  allocate(Gen_R_mo(nOrb_4times,nOrb_4times))
  allocate(Gen_H_HFB_ao(nBas_4times,nBas_4times))
  allocate(Gen_MOCoef(nBas_twice,nOrb_twice))
  allocate(Gen_MOCoef4(nBas_4times,nOrb_4times))
  allocate(Gen_U_QP(nOrb_4times,nOrb_4times))
  allocate(sw_ERI_AO(nBas_twice,nBas_twice,nBas_twice,nBas_twice))
  allocate(db_ERI_AO(nBas_twice,nBas_twice,nBas_twice,nBas_twice))

  Gen_Hc=0d0
  Gen_Hc(1:nBas           ,1:nBas           ) = Hc(1:nBas,1:nBas)
  Gen_Hc(nBas+1:nBas_twice,nBas+1:nBas_twice) = Hc(1:nBas,1:nBas)
  Gen_S=0d0
  Gen_S(1:nBas           ,1:nBas           ) = S(1:nBas,1:nBas)
  Gen_S(nBas+1:nBas_twice,nBas+1:nBas_twice) = S(1:nBas,1:nBas)
  Gen_R=0d0
  Gen_R(1:nBas                   ,1:nBas                   ) = 0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas+1:nBas_twice        ,nBas+1:nBas_twice        ) = 0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas_twice+1:nBas_3times ,nBas_twice+1:nBas_3times ) = matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb)))-0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas_3times+1:nBas_4times,nBas_3times+1:nBas_4times) = Gen_R(nBas_twice+1:nBas_3times,nBas_twice+1:nBas_3times)
  Gen_R(1:nBas                   ,nBas_3times+1:nBas_4times) =  panomMAT(1:nBas,1:nBas)
  Gen_R(nBas+1:nBas_twice        ,nBas_twice+1:nBas_3times ) = -panomMAT(1:nBas,1:nBas)
  Gen_R(nBas_twice+1:nBas_3times ,nBas+1:nBas_twice        ) = -panomMAT(1:nBas,1:nBas)
  Gen_R(nBas_3times+1:nBas_4times,1:nBas                   ) =  panomMAT(1:nBas,1:nBas)
  Gen_MOCoef=0d0
  Gen_MOCoef(1:nBas           ,1:nOrb           ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4=0d0
  Gen_MOCoef4(1:nBas                   ,1:nOrb                   ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas+1:nBas_twice        ,nOrb+1:nOrb_twice        ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas_twice+1:nBas_3times ,nOrb_twice+1:nOrb_3times ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas_3times+1:nBas_4times,nOrb_3times+1:nOrb_4times) = MOCoef(1:nBas,1:nOrb)
  sw_ERI_AO=0d0
  sw_ERI_AO(1:nBas           ,1:nBas           ,1:nBas           ,1:nBas           ) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! a a a a
  sw_ERI_AO(1:nBas           ,nBas+1:nBas_twice,1:nBas           ,nBas+1:nBas_twice) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! a b a b
  sw_ERI_AO(nBas+1:nBas_twice,1:nBas           ,nBas+1:nBas_twice,1:nBas           ) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! b a b a
  sw_ERI_AO(nBas+1:nBas_twice,nBas+1:nBas_twice,nBas+1:nBas_twice,nbas+1:nBas_twice) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! b b b b
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    do kbas=1,nBas_twice
     do lbas=1,nBas_twice
      db_ERI_AO(ibas,jbas,kbas,lbas)=sw_ERI_AO(ibas,jbas,kbas,lbas)-sw_ERI_AO(ibas,jbas,lbas,kbas)
     enddo
    enddo
   enddo
  enddo

  call Hartree_matrix_AO_basis(nBas,pMAT,ERI_AO,J)
  call exchange_matrix_AO_basis(nBas,pMAT,ERI_AO,K)
  call anomalous_matrix_AO_basis(nBas,sigma,panomMAT,ERI_AO,Delta)
  Fock(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
  Gen_H_HFB_ao(:,:)    = 0d0
  Gen_H_HFB_ao(1:nBas                   ,1:nBas                   ) = Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas+1:nBas_twice        ,nBas+1:nBas_twice        ) = Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_twice+1:nBas_3times ,nBas_twice+1:nBas_3times ) =-Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_3times+1:nBas_4times,nBas_3times+1:nBas_4times) =-Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(1:nBas                   ,nBas_3times+1:nBas_4times) = Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas+1:nBas_twice        ,nBas_twice+1:nBas_3times ) =-Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_twice+1:nBas_3times ,nBas+1:nBas_twice        ) =-Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_3times+1:nBas_4times,1:nBas                   ) = Delta(1:nBas,1:nBas)
  Gen_U_QP = matmul(transpose(Gen_MOCoef4),matmul(Gen_H_HFB_ao,Gen_MOCoef4))
  call diagonalize_matrix(nOrb_4times,Gen_U_QP,Gen_eQP_state)
  !Gen_R_mo(:,:)     = 0d0
  !do iorb=1,nOrb_twice
  ! Gen_R_mo(:,:) = Gen_R_mo(:,:) + matmul(Gen_U_QP(:,iorb:iorb),transpose(Gen_U_QP(:,iorb:iorb))) 
  !enddo
  !trace_1rdm = 0d0
  !do iorb=1,nOrb_twice
  ! trace_1rdm = trace_1rdm+Gen_R_mo(iorb,iorb) 
  !enddo
  deallocate(Gen_R_mo)
  deallocate(J)
  deallocate(K)
  deallocate(Gen_H_HFB_ao)
  deallocate(Gen_MOCoef4)
  deallocate(sw_ERI_AO)
  call scGGF2B_AO_itau_iw(nBas_twice,nBas_4times,nOrb_twice,nOrb_4times,maxSCF_GF,thresh_GF,max_diis_GF,restart_scGF2,        &
                          verbose_scGF2,chem_pot_scG,no_fock,ENuc,Gen_Hc,Gen_S,Gen_R,Gen_MOCoef,Gen_eQP_state,chem_pot,sigma, &
                          nfreqs,wcoord,wweight,Gen_U_QP,db_ERI_AO)
  deallocate(Gen_Hc)
  deallocate(Gen_S)
  deallocate(Gen_R)
  deallocate(Gen_MOCoef)
  deallocate(Gen_eQP_state)
  deallocate(Gen_U_QP)
  deallocate(db_ERI_AO)

end subroutine
