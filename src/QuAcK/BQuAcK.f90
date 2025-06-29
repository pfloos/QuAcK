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
  integer                       :: nOrb2,nBas2
  integer                       :: ibas,jbas,kbas,lbas,ifreq
  integer                       :: nO_
  integer                       :: ixyz
  integer                       :: iorb,jorb,korb,lorb

  double precision              :: chem_pot,Val
  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_qsGWB  ,end_qsGWB    ,t_qsGWB
  double precision              :: start_int    , end_int     , t_int

  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: eHFB_state(:)
  double precision,allocatable  :: U_QP(:,:)
  double precision,allocatable  :: cHFB(:,:)
  double precision,allocatable  :: PHF(:,:)
  double precision,allocatable  :: PanomHF(:,:)
  double precision,allocatable  :: FHF(:,:)
  double precision,allocatable  :: Delta(:,:)
  double precision              :: ERHF,EHFB
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)

  complex *16,allocatable       :: Chi0_ao_iw(:,:,:) 

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nBas2=nBas*nBas
  nO_=nO
  nOrb2=nOrb+nOrb

  allocate(eHF(nOrb))

  allocate(cHFB(nBas,nOrb))

  allocate(PHF(nBas,nBas))
  allocate(PanomHF(nBas,nBas))
  allocate(FHF(nBas,nBas))
  allocate(Delta(nBas,nBas))

  allocate(eHFB_state(nOrb2))
  allocate(U_QP(nOrb2,nOrb2))


  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  allocate(vMAT(nBas2,nBas2))

  allocate(Chi0_ao_iw(nfreqs,nBas2,nBas2))

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

!-----------------------------!
! Store v also as a 2D matrix !
!-----------------------------!
 
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas
      vMAT(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)=ERI_AO(ibas,jbas,kbas,lbas)
     enddo
    enddo
   enddo
  enddo

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

      call build_Xoiw_RHF_test(nBas,nBas2,nOrb,nO,cHFB,eHF,nfreqs,ntimes,wweight,wcoord,  &
                               vMAT,Chi0_ao_iw)
    endif

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call HFB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc,         &
             nBas,nOrb,nOrb2,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHFB,PHF,PanomHF,  &
             FHF,Delta,temperature,sigma,chem_pot_hf,chem_pot,restart_hfb,U_QP,eHFB_state)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

    ! Test Xo^HFB (i w) computing EcGM and EcRPA

    if(im_freqs .and. .true.) then

      call build_Xoiw_HFB_test(nBas,nBas2,nOrb,nOrb2,cHFB,eHFB_state,nfreqs,ntimes,wweight,wcoord,  &
                               U_QP,vMAT,Chi0_ao_iw)
    endif


  end if

!------------------------!
! qsGW Bogoliubov module !
!------------------------!

  if(doqsGW) then

    ! Continue with a HFB calculation
    call wall_time(start_qsGWB)
    call qsGWB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc,        &
               nBas,nOrb,nOrb2,nO_,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHFB,PHF,PanomHF,  &
               FHF,Delta,sigma,chem_pot,restart_hfb,U_QP,eHFB_state)
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
  deallocate(eHFB_state)
  deallocate(U_QP)
  deallocate(ERI_AO)
  deallocate(vMAT)
  deallocate(Chi0_ao_iw)

end subroutine
