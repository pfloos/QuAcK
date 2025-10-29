subroutine EQuAcK(working_dir,dotest,doeRHF,readFCIDUMP,nNuc,nBas,nOrb,nO,ENuc,ZNuc,rNuc,S,T,V,Hc,X,dipole_int_AO, &
                  maxSCF,max_diis,thresh,level_shift,guess_type,eweight,eforward)

! Restricted branch of Bogoliubov QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in)  :: working_dir
                                 
  logical,intent(in)             :: eforward
  logical,intent(in)             :: dotest
  logical,intent(in)             :: doeRHF
  logical,intent(in)             :: readFCIDUMP

  integer,intent(in)             :: nNuc,nBas,nOrb
  integer,intent(in)             :: nO
  double precision,intent(inout) :: ENuc

  double precision,intent(in)    :: ZNuc(nNuc),rNuc(nNuc,ncart)

  double precision,intent(inout) :: S(nBas,nBas)
  double precision,intent(inout) :: T(nBas,nBas)
  double precision,intent(inout) :: V(nBas,nBas)
  double precision,intent(inout) :: Hc(nBas,nBas)
  double precision,intent(inout) :: X(nBas,nOrb)
  double precision,intent(inout) :: dipole_int_AO(nBas,nBas,ncart)

  integer,intent(in)             :: maxSCF,max_diis
  integer,intent(in)             :: guess_type
  double precision,intent(in)    :: thresh,level_shift
  double precision,intent(in)    :: eweight

! Local variables

  logical                        :: file_exists

  integer                        :: iorb,jorb,korb,lorb
                                
  double precision               :: Val
  double precision               :: start_HF     ,end_HF       ,t_HF
  double precision               :: start_int    ,end_int      ,t_int
                                
  double precision,allocatable   :: eHF(:)
  double precision,allocatable   :: MOCoef(:,:)
  double precision,allocatable   :: pMAT(:,:)
  double precision,allocatable   :: Fock(:,:)
  double precision,allocatable   :: vMAT(:,:)
  double precision               :: EeleSD
  double precision,allocatable   :: dipole_int_MO(:,:,:)
  double precision,allocatable   :: ERI_AO(:,:,:,:)
  double precision,allocatable   :: ERI_MO(:,:,:,:)
                                
!

  write(*,*)
  write(*,*) '****************************'
  write(*,*) '* Ensemble Branch of QuAcK *'
  write(*,*) '****************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  allocate(eHF(nOrb))

  allocate(MOCoef(nBas,nOrb))

  allocate(pMAT(nBas,nBas))
  allocate(Fock(nBas,nBas))

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

!------------------------------!
! Ensemble Hartree-Fock module !
!------------------------------!

  if(doeRHF) then

    ! Run first a RHF calculation 
    call wall_time(start_HF)
    call RHF(dotest,.false.,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EeleSD,eHF,MOCoef,pMAT,Fock)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for eRHF = ',t_HF,' seconds'
    write(*,*)
    
  end if

! Memory deallocation
    
  deallocate(eHF)
  deallocate(MOCoef)
  deallocate(pMAT)
  deallocate(Fock)
  deallocate(ERI_AO)

end subroutine
