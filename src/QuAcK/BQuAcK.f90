subroutine BQuAcK(working_dir,dotest,doHFB,nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                    &
                  S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                                  &
                  guess_type,mix,temperature,sigma,chem_pot_hf,restart_hfb)

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doHFB

  logical,intent(in)            :: restart_hfb
  logical,intent(in)            :: chem_pot_hf
  integer,intent(in)            :: nNuc,nBas,nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: temperature,sigma

  double precision,intent(in)   :: ZNuc(nNuc),rNuc(nNuc,ncart)

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

! Local variables

  integer                       :: nOrb2
  integer                       :: ixyz

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_int, end_int, t_int
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO

  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: eHFB_state(:)
  double precision,allocatable  :: U_QP(:,:)
  double precision,allocatable  :: cHFB(:,:)
  double precision,allocatable  :: PHF(:,:)
  double precision,allocatable  :: PanomHF(:,:)
  double precision,allocatable  :: FHF(:,:)
  double precision,allocatable  :: Delta(:,:)
  double precision              :: ERHF,EHFB
  double precision,allocatable  :: ERI_AO(:,:,:,:)
  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: dipole_int_QP(:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: ERI_QP(:,:,:,:)

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nOrb2=nOrb+nOrb

  allocate(eHF(nOrb))

  allocate(cHFB(nBas,nOrb))

  allocate(PHF(nBas,nBas))
  allocate(PanomHF(nBas,nBas))
  allocate(FHF(nBas,nBas))
  allocate(Delta(nBas,nBas))

  allocate(eHFB_state(nOrb2))
  allocate(U_QP(nOrb2,nOrb2))


  allocate(ERI_QP(nOrb2,nOrb2,nOrb2,nOrb2))
  allocate(dipole_int_QP(nOrb2,nOrb2,ncart))

  allocate(ERI_AO(nBas,nBas,nBas,nBas))
  call wall_time(start_int)
  call read_2e_integrals(working_dir,nBas,ERI_AO)
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
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHFB,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call HFB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc,       &
             nBas,nOrb,nOrb2,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHFB,PHF,PanomHF,  &
             FHF,Delta,temperature,sigma,chem_pot_hf,restart_hfb,U_QP,eHFB_state)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

  call wall_time(start_AOtoMO)
  
  write(*,*)
  write(*,*) 'AO to MO transformation... Please be patient'
  write(*,*)

  ! Read and transform dipole-related integrals
  
  allocate(dipole_int_MO(nOrb,nOrb,ncart))
  do ixyz=1,ncart
    call AOtoMO(nBas,nOrb,cHFB,dipole_int_AO(1,1,ixyz),dipole_int_MO(1,1,ixyz))
    call AOtoMO(nOrb,nOrb2,U_QP,dipole_int_MO(1,1,ixyz),dipole_int_QP(1,1,ixyz))
  end do 
  deallocate(dipole_int_MO)
  
  ! 4-index transform 
  
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,cHFB,ERI_AO,ERI_MO)
  call AOtoMO_ERI_RHF(nOrb,nOrb2,U_QP,ERI_MO,ERI_QP)
  deallocate(ERI_MO)

  call wall_time(end_AOtoMO)
  
  t_AOtoMO = end_AOtoMO - start_AOtoMO
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for AO to MO transformation = ',t_AOtoMO,' seconds'
  write(*,*)


! Memory deallocation
    
  deallocate(eHF)
  deallocate(cHFB)
  deallocate(PHF)
  deallocate(PanomHF)
  deallocate(FHF)
  deallocate(Delta)
  deallocate(ERI_AO)
  deallocate(ERI_QP)
  deallocate(dipole_int_QP)
  deallocate(eHFB_state)
  deallocate(U_QP)

end subroutine
