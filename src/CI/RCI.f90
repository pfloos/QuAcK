
! ---

subroutine RCI(dotest, doCIS, doCIS_D, doCID, doCISD, doFCI, singlet, triplet, nOrb, &
               nC, nO, nV, nR, nS, ERI, dipole_int, epsHF, EHF)

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
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  double precision              :: start_CI     ,end_CI       ,t_CI

!------------------------------------------------------------------------
! Compute CIS excitations
!------------------------------------------------------------------------

  if(doCIS) then

    call wall_time(start_CI)
    call RCIS(dotest,singlet,triplet,doCIS_D,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,epsHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CIS = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CID excitations
!------------------------------------------------------------------------

  if(doCID) then

    call wall_time(start_CI)
    call CID(dotest,singlet,triplet,nOrb,nC,nO,nV,nR,ERI,epsHF,EHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CID = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CISD excitations
!------------------------------------------------------------------------

  if(doCISD) then

    call wall_time(start_CI)
    call CISD(dotest,singlet,triplet,nOrb,nC,nO,nV,nR,ERI,epsHF,EHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CISD = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute FCI 
!------------------------------------------------------------------------

  if(doFCI) then

    call wall_time(start_CI)
    write(*,*) ' FCI is not yet implemented! Sorry.'
!   call FCI(nOrb,nC,nO,nV,nR,ERI,epsHF)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for FCI = ',t_CI,' seconds'
    write(*,*)

  end if

end subroutine
