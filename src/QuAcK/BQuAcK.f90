subroutine BQuAcK(working_dir,dotest,doHFB,nNuc,nBas,nOrb,nC,nO,nV,nR,ENuc,ZNuc,rNuc,                                    &
                  S,T,V,Hc,X,dipole_int_AO,maxSCF_HF,max_diis_HF,thresh_HF,level_shift,                                  &
                  guess_type,mix)

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doHFB

  integer,intent(in)            :: nNuc,nBas,nOrb
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
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)            :: maxSCF_HF,max_diis_HF
  double precision,intent(in)   :: thresh_HF,level_shift,mix
  integer,intent(in)            :: guess_type

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF

  double precision              :: start_int, end_int, t_int
  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: cHF(:,:)
  double precision,allocatable  :: PHF(:,:)
  double precision,allocatable  :: PanomHF(:,:)
  double precision,allocatable  :: FHF(:,:)
  double precision              :: ERHF,EHFB
!  double precision,allocatable  :: dipole_int_MO(:,:,:)
  double precision,allocatable  :: ERI_AO(:,:,:,:)
!  double precision,allocatable  :: ERI_MO(:,:,:,:)

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  allocate(eHF(nOrb))
  allocate(cHF(nBas,nOrb))
  allocate(PHF(nBas,nBas))
  allocate(PanomHF(nBas,nBas))
  allocate(FHF(nBas,nBas))
!  allocate(dipole_int_MO(nOrb,nOrb,ncart))
!  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))

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
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,eHF,cHF,PHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

    ! Continue with a HFB calculation
    call wall_time(start_HF)
    call HFB(dotest,maxSCF_HF,thresh_HF,max_diis_HF,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHFB,eHF,cHF,PHF,PanomHF,FHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

! Memory deallocation
    
  deallocate(eHF)
  deallocate(cHF)
  deallocate(PHF)
  deallocate(PanomHF)
  deallocate(FHF)
!  deallocate(dipole_int_MO)
!  deallocate(ERI_MO)
  deallocate(ERI_AO)

end subroutine
