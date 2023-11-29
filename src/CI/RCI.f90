subroutine RCI(dotest,doCIS,doCIS_D,doCID,doCISD,doFCI,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,dipole_int, &
               epsHF,EHF,cHF,S)

! Configuration interaction module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doCIS
  logical,intent(in)            :: doCIS_D
  logical,intent(in)            :: doCID
  logical,intent(in)            :: doCISD
  logical,intent(in)            :: doFCI

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
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
    call RCIS(dotest,singlet,triplet,doCIS_D,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,epsHF)
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
    call CID(dotest,singlet,triplet,nBas,nC,nO,nV,nR,ERI,epsHF,EHF)
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
    call CISD(dotest,singlet,triplet,nBas,nC,nO,nV,nR,ERI,epsHF,EHF)
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
