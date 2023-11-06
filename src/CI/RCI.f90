subroutine RCI(doCIS,doCIS_D,doCID,doCISD,doFCI,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,dipole_int, &
               epsHF,EHF,cHF,S)

! Configuration interaction module

  implicit none
  include 'parameters.h'

! Input variables

  logical                       :: doCIS
  logical                       :: doCIS_D
  logical                       :: doCID
  logical                       :: doCISD
  logical                       :: doFCI

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_CI     ,end_CI       ,t_CI

!------------------------------------------------------------------------
! Compute CIS excitations
!------------------------------------------------------------------------

  if(doCIS) then

    call wall_time(start_CI)
    call CIS(singlet,triplet,doCIS_D,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,epsHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CID excitations
!------------------------------------------------------------------------

  if(doCID) then

    call wall_time(start_CI)
    call CID(singlet,triplet,nBas,nC,nO,nV,nR,ERI,epsHF,EHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CID = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CISD excitations
!------------------------------------------------------------------------

  if(doCISD) then

    call wall_time(start_CI)
    call CISD(singlet,triplet,nBas,nC,nO,nV,nR,ERI,epsHF,EHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CISD = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute FCI 
!------------------------------------------------------------------------

  if(doFCI) then

    call wall_time(start_CI)
    write(*,*) ' FCI is not yet implemented! Sorry.'
!   call FCI(nBas,nC,nO,nV,nR,ERI,epsHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for FCI = ',t_CI,' seconds'
    write(*,*)

  end if

end subroutine