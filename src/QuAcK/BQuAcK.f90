subroutine BQuAcK(working_dir,dotest,doHFB,dophRPA,doqsGW,nNuc,nBas,nOrb,nO,ENuc,eta,shift,ZNuc,   &
                  rNuc,S,T,V,Hc,X,dipole_int_AO,maxSCF,max_diis,thresh,level_shift,guess_type,mix, &
                  temperature,sigma,chem_pot_hf,restart_hfb,nfreqs,ntimes,wcoord,wweight)

! Restricted branch of Bogoliubov QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in)  :: working_dir
                                 
  logical,intent(in)             :: dotest
                                 
  logical,intent(in)             :: doHFB
  logical,intent(in)             :: dophRPA
  logical,intent(in)             :: doqsGW

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
  double precision,intent(in)    :: thresh,level_shift,mix
  integer,intent(in)             :: guess_type

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
  double precision,allocatable   :: panomMAT(:,:)
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
  endif
  close(314)
  
  call wall_time(end_int)
  t_int = end_int - start_int

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for reading 2e-integrals =',t_int,' seconds'
  write(*,*)

!--------------------------------!
! Hartree-Fock Bogoliubov module !
!--------------------------------!

  if(doHFB) then

    ! Run first a RHF calculation 
    call wall_time(start_HF)
    call RHF(dotest,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)
    
    ! An imaginary frequencies implementation of qsGW for testing (expensive!)
    inquire(file='qsRGWim', exist=file_exists)
    if(doqsGW .and. file_exists) then

      ! Continue with a HFB calculation
      call wall_time(start_qsGWB)
      call qsRGWim(dotest,maxSCF,thresh,max_diis,level_shift,eta,shift,nNuc,ZNuc,rNuc,ENuc, &
                   nBas,nOrb,nO_,verbose,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef, &
                   pMAT,Fock,nfreqs,ntimes,wcoord,wweight)
      call wall_time(end_qsGWB)

      t_qsGWB = end_qsGWB - start_qsGWB
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGWB = ',t_qsGWB,' seconds'
      write(*,*)

    end if

    ! Compute EcRPA and EcGM energies for RHF
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
     ! Test down-folded G0W0 matrix?
     !call dfRG0W0mat(nOrb,nO,eta,shift,eHF,vMAT,nfreqs,ntimes,wcoord,wweight)
     ! Test EcGM computed from Sigma_c(iw)
     !call EcGM_w_RHF_Sigma(nOrb,nO,1,eHF,nfreqs,wweight,wcoord,vMAT,EeleSD+Enuc,EcGM)
     deallocate(vMAT)
     call wall_time(end_Ecorr)

     t_Ecorr = end_Ecorr - start_Ecorr
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
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------!
! qsGW Bogoliubov module !
!------------------------!

  if(doqsGW) then

    ! Continue with a HFB calculation
    call wall_time(start_qsGWB)
    call qsRGWBim(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,eta,shift,    &
               nBas,nOrb,nOrb_twice,nO_,verbose,S,T,V,Hc,ERI_AO,dipole_int_AO,X,Eelec,        & 
               MOCoef,pMAT,panomMAT,Fock,Delta,sigma,chem_pot,U_QP,eQP_state,nfreqs,ntimes,   &
               wcoord,wweight) 
    call wall_time(end_qsGWB)

    t_qsGWB = end_qsGWB - start_qsGWB
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGWB = ',t_qsGWB,' seconds'
    write(*,*)

  end if

  ! Compute EcRPA and EcGM energies for qsGWB or HFB
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
   call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                          U_QP,Eelec+ENuc,EcRPA,EcGM)
   ! Test down-folded G0W0 Bogoliubov matrix?
   !call dfRG0W0Bmat(nOrb,nOrb_twice,eta,shift,eQP_state,U_QP,vMAT,nfreqs,ntimes,wcoord,wweight)
   ! Test EcGM computed from Sigma_c(iw)
   !call EcGM_w_RHFB_Sigma(nOrb,nOrb_twice,1,eQP_state,nfreqs,wweight,wcoord,vMAT,U_QP,EeleSD+Enuc,EcGM)
   deallocate(vMAT)
   call wall_time(end_Ecorr)

   t_Ecorr = end_Ecorr - start_Ecorr
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
