subroutine BQuAcK(working_dir,dotest,doHFB,doqsGW,nNuc,nBas,nOrb,nO,ENuc,ZNuc,rNuc,                   &
                  S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,               &
                  guess_type,mix,temperature,sigma,chem_pot_hf,restart_hfb,im_freqs,nfreqs,ntimes,wcoord,wweight)

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doHFB
  logical,intent(in)            :: doqsGW

  logical,intent(in)            :: restart_hfb
  logical,intent(in)            :: im_freqs
  logical,intent(in)            :: chem_pot_hf
  integer,intent(in)            :: nNuc,nBas,nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  double precision,intent(inout):: ENuc
  double precision,intent(in)   :: temperature,sigma

  double precision,intent(in)   :: ZNuc(nNuc),rNuc(nNuc,ncart)
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)

  double precision,intent(inout)   :: S(nBas,nBas)
  double precision,intent(inout)   :: T(nBas,nBas)
  double precision,intent(inout)   :: V(nBas,nBas)
  double precision,intent(inout)   :: Hc(nBas,nBas)
  double precision,intent(inout)   :: X(nBas,nOrb)
  double precision,intent(inout)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

! Local variables

  logical                       :: file_exists
  integer                       :: nOrb_twice
  integer                       :: nO_
  integer                       :: ixyz
  integer                       :: iorb,jorb,korb,lorb

  double precision              :: chem_pot,Val
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_qsGWB  ,end_qsGWB    ,t_qsGWB
  double precision              :: start_int    , end_int     , t_int

  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: eONEBODY_state(:)
  double precision,allocatable  :: U_QP(:,:)
  double precision,allocatable  :: cHFB(:,:)
  double precision,allocatable  :: PHF(:,:)
  double precision,allocatable  :: PanomHF(:,:)
  double precision,allocatable  :: FHF(:,:)
  double precision,allocatable  :: Delta(:,:)
  double precision              :: ERHF,EHFB
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)

!

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nO_=nO
  nOrb_twice=nOrb+nOrb

  allocate(eHF(nOrb))

  allocate(cHFB(nBas,nOrb))

  allocate(PHF(nBas,nBas))
  allocate(PanomHF(nBas,nBas))
  allocate(FHF(nBas,nBas))
  allocate(Delta(nBas,nBas))

  allocate(eONEBODY_state(nOrb_twice))
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
    call RHF(dotest,maxSCF_HF,thresh_HF,max_diis_HF,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHFB,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

    ! Test Xo^HF (i w) computing EcGM and EcRPA
    if(im_freqs .and. .true.) then

      call Xoiw_RHF_tests(nBas,nOrb,nO,cHFB,eHF,nfreqs,ntimes,wweight,wcoord,ERI_AO)

    endif

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call HFB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc,              &
             nBas,nOrb,nOrb_twice,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHFB,PHF,PanomHF,  &
             FHF,Delta,temperature,sigma,chem_pot_hf,chem_pot,restart_hfb,U_QP,eONEBODY_state)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

    ! Test Xo^HFB (i w) computing EcGM and EcRPA
    if(im_freqs .and. .true.) then

      call Xoiw_HFB_tests(nBas,nOrb,nOrb_twice,cHFB,eONEBODY_state,nfreqs,ntimes,wweight,wcoord,U_QP,ERI_AO)

    endif

  end if

!------------------------!
! qsGW Bogoliubov module !
!------------------------!

  if(doqsGW) then

    ! Continue with a HFB calculation
    call wall_time(start_qsGWB)
    call qsGWB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc,              &
               nBas,nOrb,nOrb_twice,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHFB,PHF,PanomHF,  &
               FHF,Delta,sigma,chem_pot,restart_hfb,U_QP,eONEBODY_state)
    call wall_time(end_qsGWB)

    t_qsGWB = end_qsGWB - start_qsGWB
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGWB = ',t_qsGWB,' seconds'
    write(*,*)

  end if

! Memory deallocation
    
  deallocate(eHF)
  deallocate(cHFB)
  deallocate(PHF)
  deallocate(PanomHF)
  deallocate(FHF)
  deallocate(Delta)
  deallocate(eONEBODY_state)
  deallocate(U_QP)
  deallocate(ERI_AO)

end subroutine
