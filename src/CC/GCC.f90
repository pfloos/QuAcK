subroutine GCC(dotest,doCCD,dopCCD,doDCD,doCCSD,doCCSDT,dodrCCD,dorCCD,docrCCD,dolCCD, & 
               maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)

! Generalized Coupled-cluster module

  implicit none
  include 'parameters.h'

! Input variables

  logical                       :: dotest

  logical                       :: doCCD
  logical                       :: dopCCD
  logical                       :: doDCD
  logical                       :: doCCSD
  logical                       :: doCCSDT
  logical                       :: dodrCCD
  logical                       :: dorCCD
  logical                       :: docrCCD
  logical                       :: dolCCD

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: start_CC     ,end_CC       ,t_CC

!------------------------------------------------------------------------
! Perform CCD calculation
!------------------------------------------------------------------------

  if(doCCD) then

    call wall_time(start_CC)
    call GCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform DCD calculation
!------------------------------------------------------------------------

  if(doDCD) then

    call wall_time(start_CC)
!   call DCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR, & 
!            ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for DCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform CCSD or CCSD(T) calculation
!------------------------------------------------------------------------

  if(doCCSDT) doCCSD = .true.

  if(doCCSD) then

    call wall_time(start_CC)
    call GCCSD(dotest,maxSCF,thresh,max_diis,doCCSDT,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCSD or CCSD(T)= ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform direct ring CCD calculation
!------------------------------------------------------------------------

  if(dodrCCD) then

    call wall_time(start_CC)
    call drGCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for direct ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ring CCD calculation
!------------------------------------------------------------------------

  if(dorCCD) then

    call wall_time(start_CC)
    call rGCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for rCCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform crossed-ring CCD calculation
!------------------------------------------------------------------------

  if(docrCCD) then

    call wall_time(start_CC)
    call crGCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for crossed-ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ladder CCD calculation
!------------------------------------------------------------------------

  if(dolCCD) then

    call wall_time(start_CC)
    call lGCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ladder CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform pair CCD calculation
!------------------------------------------------------------------------

  if(dopCCD) then

    call wall_time(start_CC)
!   call pCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pair CCD = ',t_CC,' seconds'
    write(*,*)

  end if

end subroutine
