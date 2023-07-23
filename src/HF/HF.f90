subroutine HF(doRHF,doUHF,doRMOM,doUMOM,unrestricted,maxSCF,thresh,max_diis,guess_type,mix,level_shift, & 
              nNuc,ZNuc,rNuc,ENuc,nBas,nO,S,T,V,Hc,F,ERI,dipole_int,X,EHF,epsHF,cHF,PHF)

! Hartree-Fock module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doRHF
  logical,intent(in)            :: doUHF
  logical,intent(in)            :: doRMOM
  logical,intent(in)            :: doUMOM

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift
  logical,intent(in)            :: mix

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF

  integer                       :: nSCF
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: Gap 
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error(:,:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)

! Output variables

  logical,intent(out)           :: unrestricted

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: epsHF(nBas)
  double precision,intent(out)  :: cHF(nBas,nBas)
  double precision,intent(out)  :: PHF(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

!------------------------------------------------------------------------
! Compute RHF energy
!------------------------------------------------------------------------

  if(doRHF) then

   ! Check that RHF calculation is worth doing...

    if(nO(1) /= nO(2)) then
      write(*,*) ' !!! The system does not appear to be closed shell !!!'
      write(*,*)
      stop
    end if

    call wall_time(start_HF)
    call RHF(maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,F,ERI,dipole_int,X,EHF,epsHF,cHF,PHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HF = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute UHF energy
!------------------------------------------------------------------------

  if(doUHF) then

    ! Switch on the unrestricted flag
    unrestricted = .true.

    call wall_time(start_HF)
    call UHF(maxSCF,thresh,max_diis,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI,dipole_int,X,EHF,epsHF,cHF,PHF)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HF = ',t_HF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Restricted maximum overlap method
!------------------------------------------------------------------------

  if(doRMOM) then

    ! Check that RMOM calculation is worth doing...

    if(nO(1) /= nO(2)) then
      write(*,*) ' !!! The system does not appear to be closed shell !!!'
      write(*,*)
      stop
    end if

!   call RMOM(maxSCF,thresh,max_diis,guess_type,nNuc,ZNuc,rNuc,ENuc, &
!           nBas,nO,S,T,V,Hc,ERI,dipole_int,X,EHF,epsHF,cHF,PHF)

  end if

!------------------------------------------------------------------------
! Unrestricted maximum overlap method
!------------------------------------------------------------------------

  if(doUMOM) then

    ! Switch on the unrestricted flag
    unrestricted = .true.

!   call UMOM(maxSCF,thresh,max_diis,guess_type,nNuc,ZNuc,rNuc,ENuc, &
!           nBas,nO,S,T,V,Hc,ERI,dipole_int,X,EHF,epsHF,cHF,PHF)

  end if

end subroutine 
