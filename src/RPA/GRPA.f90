subroutine GRPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

! Random-phase approximation module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: dophRPA
  logical,intent(in)            :: dophRPAx
  logical,intent(in)            :: docrRPA
  logical,intent(in)            :: doppRPA

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(dophRPA) then

    call wall_time(start_RPA)
    call phGRPA(dotest,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute RPAx (RPA with exchange) excitations
!------------------------------------------------------------------------

  if(dophRPAx) then

    call wall_time(start_RPA)
    call phGRPAx(dotest,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPAx = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute crRPA excitations
!------------------------------------------------------------------------

  if(docrRPA) then

    call wall_time(start_RPA)
    call crGRPA(dotest,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for cr-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ppRPA excitations
!------------------------------------------------------------------------

  if(doppRPA) then

    call wall_time(start_RPA)
    call ppGRPA(dotest,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

end subroutine
