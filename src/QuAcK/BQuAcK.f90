subroutine BQuAcK(working_dir,dotest,doRHFB,doBRPA,dophRPA,doG0W0,doqsGW,nNuc,nBas,nOrb,nO,ENuc,eta,shift, &
                  ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO,maxSCF,max_diis,thresh,level_shift,guess_type,        &
                  maxSCF_GW,max_diis_GW,thresh_GW,dolinGW,                                                 &
                  temperature,sigma,chem_pot_hf,restart_hfb,nfreqs,ntimes,wcoord,wweight)

! Restricted branch of Bogoliubov QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in)  :: working_dir
                                 
  logical,intent(in)             :: dotest
                                 
  logical,intent(in)             :: doRHFB
  logical,intent(in)             :: doBRPA
  logical,intent(in)             :: dophRPA
  logical,intent(in)             :: doG0W0
  logical,intent(in)             :: doqsGW
  logical,intent(in)             :: dolinGW

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

  logical                        :: file_exists
  integer                        :: verbose
  integer                        :: nOrb_twice
  integer                        :: nO_
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
  nO_=nO
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

! For the Hubbard model read two-body parameters (U, J, etc from hubbard file)

  inquire(file='hubbard', exist=file_exists)
  if(file_exists) then
   write(*,*)
   write(*,*) 'Reading Hubbard model two-body parameters'
   write(*,*)
   ERI_AO=0d0; ENuc=0d0;
   open(unit=314, form='formatted', file='hubbard', status='old')
   do
    read(314,*) iorb,jorb,korb,lorb,Val
    if(iorb==jorb .and. jorb==korb .and. korb==lorb .and. iorb==-1) then
     nO_ = int(Val)
     cycle
    endif
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
    call RHF(dotest,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)
    
    ! An imaginary frequencies implementation of qsGW for testing (expensive!)
    inquire(file='qsRGWim', exist=file_exists)
    if(doqsGW .and. file_exists) then

      ! Continue with a HFB calculation
      call wall_time(start_qsGWB)
      call qsRGWim(dotest,maxSCF_GW,thresh_GW,max_diis_GW,level_shift,eta,shift,nNuc,ZNuc,rNuc,ENuc, &
                   nBas,nOrb,nO_,verbose,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,          &
                   pMAT,Fock,nfreqs,ntimes,wcoord,wweight)
      call wall_time(end_qsGWB)

      t_qsGWB = end_qsGWB - start_qsGWB
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGWB = ',t_qsGWB,' seconds'
      write(*,*)

    end if

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
     ! Test down-folded G0W0 matrix?
     if(doG0W0) then
      call dfRG0W0mat(nOrb,nO,eta,shift,eHF,vMAT,nfreqs,ntimes,wcoord,wweight)
     endif
     ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr
     if(dolinGW .and. dophRPA) then
      allocate(pMATcorr(nBas,nBas))
      call linDyson_GW_RHF(nBas,nOrb,nO,MOCoef,eHF,nfreqs,wweight,wcoord,ERI_AO,vMAT, &
                          Enuc,EcGM,T,V,S,pMAT,pMATcorr)
      !call quadDyson_GW_RHF(nBas,nOrb,nO,MOCoef,eHF,nfreqs,wweight,wcoord,ERI_AO,vMAT, &
      !                    Enuc,EcGM,T,V,S,pMAT,pMATcorr)
      !call cubDyson_GW_RHF(nBas,nOrb,nO,MOCoef,eHF,nfreqs,wweight,wcoord,ERI_AO,vMAT, &
      !                    Enuc,EcGM,T,V,S,pMAT,pMATcorr)
      deallocate(pMATcorr)
     endif
     ! Test EcGM computed from Sigma_c(iw) [ NOTE: This is really bad numerically and never used. ]
     !call EcGM_w_RHF_Sigma(nOrb,nO,1,eHF,nfreqs,wweight,wcoord,vMAT,EeleSD+Enuc,EcGM)
     deallocate(vMAT)
     call wall_time(end_Ecorr)

     t_Ecorr = end_Ecorr - start_Ecorr
     write(*,*)
     write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Ecorr = ',t_Ecorr,' seconds'
     write(*,*)

    endif

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call RHFB(dotest,doqsGW,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,       &
             nBas,nOrb,nOrb_twice,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,Eelec,eHF,MOCoef,   & 
             pMAT,panomMAT,Fock,Delta,temperature,sigma,chem_pot_hf,chem_pot,restart_hfb, &
             U_QP,eQP_state)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------!
! qsGW Bogoliubov module !
!------------------------!

  if(doqsGW) then

    ! Continue with a HFB calculation
    call wall_time(start_qsGWB)
    call qsRGWBim(dotest,maxSCF_GW,thresh_GW,max_diis_GW,level_shift,nNuc,ZNuc,rNuc,ENuc,eta,shift, &
               nBas,nOrb,nOrb_twice,nO_,verbose,S,T,V,Hc,ERI_AO,dipole_int_AO,X,Eelec,              & 
               MOCoef,pMAT,panomMAT,Fock,Delta,sigma,chem_pot,U_QP,eQP_state,nfreqs,ntimes,         &
               wcoord,wweight) 
    call wall_time(end_qsGWB)

    t_qsGWB = end_qsGWB - start_qsGWB
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGWB = ',t_qsGWB,' seconds'
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
   ! Test down-folded G0W0 Bogoliubov matrix?
   if(doG0W0) then
    call dfRG0W0Bmat(nBas,nO_,nOrb,nOrb_twice,chem_pot,eta,shift,eQP_state,U_QP,vMAT,nfreqs,ntimes, &
                     wcoord,wweight,sigma,S,T,V,Hc,MOCoef,pMAT,panomMAT,Delta,ERI_AO)
   endif
   ! Test linearized-Dyson equation G ~ Go + Go Sigma_c Go -> Pcorr and Panomcorr
   if(dolinGW .and. dophRPA) then
    allocate(pMATcorr(nBas,nBas),panomMATcorr(nBas,nBas))
    call linDyson_GW_RHFB(nBas,nOrb,nOrb_twice,MOCoef,eQP_state,nfreqs,wweight,wcoord,ERI_AO,vMAT,U_QP,&
                          Enuc,EcGM,sigma,T,V,S,pMAT,panomMAT,pMATcorr,panomMATcorr)
    !call quadDyson_GW_RHFB(nBas,nOrb,nOrb_twice,MOCoef,eQP_state,nfreqs,wweight,wcoord,ERI_AO,vMAT,U_QP,&
    !                      Enuc,EcGM,sigma,T,V,S,pMAT,panomMAT,pMATcorr,panomMATcorr)
    !call cubDyson_GW_RHFB(nBas,nOrb,nOrb_twice,MOCoef,eQP_state,nfreqs,wweight,wcoord,ERI_AO,vMAT,U_QP,&
    !                      Enuc,EcGM,sigma,T,V,S,pMAT,panomMAT,pMATcorr,panomMATcorr)
    deallocate(pMATcorr,panomMATcorr)
   endif
   ! Test EcGM computed from Sigma_c(iw) [ NOTE: This is really bad numerically and never used. ]
   !call EcGM_w_RHFB_Sigma(nOrb,nOrb_twice,1,eQP_state,nfreqs,wweight,wcoord,vMAT,U_QP,EeleSD+Enuc,EcGM)
   deallocate(vMAT)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
   write(*,*)
   write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Ecorr = ',t_Ecorr,' seconds'
   write(*,*)

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
   call EcRPA_Bogoliubov_FCIDUMP(nO_,nOrb,nOrb_twice,sigma,ERI_MO,vMAT,Fock,Delta,pMAT,panomMAT,eQP_state,U_QP, &
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

subroutine EcRPA_Bogoliubov_FCIDUMP(nO,nOrb,nOrb_twice,sigma,ERI_MO,vMAT,Fock,Delta,pMAT,panomMAT,eQP_state,U_QP, &
                                    chem_pot,ntimes,nfreqs,wcoord,wweight)
  implicit none
  include 'parameters.h'
!
  integer,intent(in)             :: nO
  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nOrb_twice
  integer,intent(in)             :: ntimes
  integer,intent(in)             :: nfreqs
  double precision,intent(in)    :: sigma
  double precision,intent(in)    :: wcoord(nfreqs),wweight(nfreqs)

! Local variables

  logical                        :: file_exists
  integer                        :: iorb,jorb,korb,lorb
  double precision               :: Ecore,ENuc,EJ,EK,EL,Eelec,EHFB,thrs_N,trace_1rdm,Val,EcRPA,EcGM  
  double precision,external      :: trace_matrix
  double precision,allocatable   :: J(:,:)
  double precision,allocatable   :: K(:,:)
  double precision,allocatable   :: H_HFB(:,:)
  double precision,allocatable   :: Hc(:,:)
  double precision,allocatable   :: R(:,:)

! Output variables

  double precision,intent(inout) :: chem_pot
  double precision,intent(out)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(out)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(out)   :: Fock(nOrb,nOrb)
  double precision,intent(out)   :: Delta(nOrb,nOrb)
  double precision,intent(out)   :: pMAT(nOrb,nOrb)
  double precision,intent(out)   :: panomMAT(nOrb,nOrb)
  double precision,intent(out)   :: eQP_state(nOrb_twice)
  double precision,intent(out)   :: U_QP(nOrb_twice,nOrb_twice)

  character(len=100)             :: arg                              
!

  write(*,*)
  write(*,*) '*********************************'
  write(*,*) '* Bogoliubov EcRPA from FCIDUMP *'
  write(*,*) '*********************************'
  write(*,*)

  ENuc=0d0
  allocate(J(nOrb,nOrb))
  allocate(K(nOrb,nOrb))
  allocate(Hc(nOrb,nOrb))
  allocate(H_HFB(nOrb_twice,nOrb_twice))
  allocate(R(nOrb_twice,nOrb_twice))
  
  inquire(file='FCIDUMP', exist=file_exists)
  if(file_exists) then
   write(*,*) 'Reading FCIDUMP file before computing EcRPA Bogoliubov'
   ERI_MO=0d0; vMAT=0d0; ENuc=0d0; Hc=0d0;
   open(unit=314, form='formatted', file='FCIDUMP', status='old')
   read(314,'(a)') arg
   read(314,'(a)') arg
   read(314,'(a)') arg
   do
    read(314,*) Val,iorb,jorb,korb,lorb
    if(iorb==jorb .and. jorb==korb .and. korb==lorb .and. iorb==0) then
     ENuc=Val
     exit
    endif
    if(korb==lorb .and. lorb==0) then
     Hc(iorb,jorb)=Val
     Hc(jorb,iorb)=Hc(iorb,jorb)
    else
     ERI_MO(iorb,korb,jorb,lorb)=Val
     ERI_MO(iorb,lorb,jorb,korb)=Val
     ERI_MO(jorb,lorb,iorb,korb)=Val
     ERI_MO(jorb,korb,iorb,lorb)=Val
     ERI_MO(korb,iorb,lorb,jorb)=Val
     ERI_MO(lorb,iorb,korb,jorb)=Val
     ERI_MO(lorb,jorb,korb,iorb)=Val
     ERI_MO(korb,jorb,lorb,iorb)=Val
    endif
   enddo
   close(314)
  endif
  do iorb=1,nOrb
   do jorb=1,nOrb
    do korb=1,nOrb
     do lorb=1,nOrb
      vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
     enddo
    enddo
   enddo
  enddo

   inquire(file='occupancies', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading occupancies file before computing EcRPA Bogoliubov'
    pMAT=0d0
    panomMAT=0d0
    open(unit=314, form='formatted', file='occupancies', status='old')
    do iorb=1,nOrb
      read(314,*) Val
      pMAT(iorb,iorb)=Val
      Val=0.5d0*Val
      panomMAT(iorb,iorb)=sqrt(abs(Val*(1.0d0-Val)))
    enddo
    close(314)
   endif
   
   call Hartree_matrix_AO_basis(nOrb,pMAT,ERI_MO,J)
   call exchange_matrix_AO_basis(nOrb,pMAT,ERI_MO,K)
   call anomalous_matrix_AO_basis(nOrb,sigma,panomMAT,ERI_MO,Delta)
   thrs_N=1d-6
   Fock(:,:)=Hc(:,:)+J(:,:)+0.5d0*K(:,:)
   do iorb=1,nOrb
   Fock(iorb,iorb)=Fock(iorb,iorb)-chem_pot
   enddo
   ! Hcore energy
   Ecore = trace_matrix(nOrb,matmul(pMAT,Hc))
   ! Coulomb energy
   EJ = 0.5d0*trace_matrix(nOrb,matmul(pMAT,J))
   ! Exchange energy
   EK = 0.25d0*trace_matrix(nOrb,matmul(pMAT,K))
   ! Anomalous energy
   EL = trace_matrix(nOrb,matmul(panomMAT,Delta))
   ! Total electronic energy
   Eelec = Ecore + EJ + EK + EL
   ! Total energy
   EHFB = Eelec + ENuc
   H_HFB(1:nOrb,1:nOrb)=Fock(1:nOrb,1:nOrb)
   H_HFB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)=-Fock(1:nOrb,1:nOrb)
   H_HFB(1:nOrb,nOrb+1:nOrb_twice)=Delta(1:nOrb,1:nOrb)
   H_HFB(nOrb+1:nOrb_twice,1:nOrb)=Delta(1:nOrb,1:nOrb)
   call diagonalize_matrix(nOrb_twice,U_QP,eQP_state)
   ! Build R 
   R(:,:)     = 0d0
   do iorb=1,nOrb
    R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb)))
   enddo
   trace_1rdm = 0d0
   do iorb=1,nOrb
    trace_1rdm = trace_1rdm+R(iorb,iorb)
   enddo
   ! Adjust the chemical potential 
   if( abs(trace_1rdm-nO) > thrs_N ) &
    call fix_chem_pot(nO,nOrb,nOrb_twice,0,thrs_N,trace_1rdm,chem_pot,H_HFB,U_QP,R,eQP_state)
 
   write(*,*)
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33)')           ' Summary               '
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',Ecore,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK + EL,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Anomalous    energy = ',EL,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',Eelec,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
   write(*,'(A33,1X,F16.10,A3)') ' HFB          energy = ',EHFB,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,*)
   write(*,'(A50)') '---------------------------------------'
   write(*,'(A50)') ' HFB QP energies '
   write(*,'(A50)') '---------------------------------------'
   do iorb=1,nOrb_twice
    write(*,'(I7,10F15.8)') iorb,eQP_state(iorb)
   enddo
   write(*,*)
   
   !U_QP(:,1:nOrb)=0d0
   !do iorb=1,nOrb
   ! U_QP(iorb,iorb) = sqrt(abs(0.5d0*pMAT(iorb,iorb)))
   ! U_QP(iorb+nOrb,iorb) = sqrt(abs(1.0d0-0.5d0*pMAT(iorb,iorb)))
   !enddo
   Eelec = 0d0
   call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                          U_QP,EHFB,EcRPA,EcGM)
   deallocate(J,K,Hc)
   deallocate(H_HFB,R)

end subroutine

