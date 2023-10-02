subroutine CC(doCCD,dopCCD,doDCD,doCCSD,doCCSDT,dodrCCD,dorCCD,docrCCD,dolCCD, & 
              maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)

! Coupled-cluster module

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: start_CC     ,end_CC       ,t_CC

!------------------------------------------------------------------------
! Perform CCD calculation
!------------------------------------------------------------------------

  if(doCCD) then

    call wall_time(start_CC)
    call CCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
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
    call DCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR, & 
             ERI,ENuc,EHF,epsHF)
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
    call CCSD(maxSCF,thresh,max_diis,doCCSDT,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
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
    call drCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
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
    call rCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF,epsHF)
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
    call crCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
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
    call lCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
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
!   call pCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call ROpCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)

    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pair CCD = ',t_CC,' seconds'
    write(*,*)

  end if

end subroutine
