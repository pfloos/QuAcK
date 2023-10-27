subroutine GRPA(dophRPA,dophRPAx,doppRPA,TDA,doACFDT,exchange_kernel,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)

! Random-phase approximation module

  implicit none
  include 'parameters.h'

! Input variables

  logical                       :: dophRPA
  logical                       :: dophRPAx
  logical                       :: doppRPA

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(dophRPA) then

    call wall_time(start_RPA)
    call phGRPA(TDA,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)
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
    call phGRPAx(TDA,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)

    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPAx = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ppRPA excitations
!------------------------------------------------------------------------

  if(doppRPA) then

    call wall_time(start_RPA)
    call ppGRPA(TDA,doACFDT,nBas,nC,nO,nV,nR,ENuc,EHF,ERI,dipole_int,epsHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

end subroutine
