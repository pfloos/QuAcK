
! ---

subroutine RCC(dotest, doCCD, dopCCD, doDCD, doCCSD, doCCSDT, dodrCCD, dorCCD, docrCCD, dolCCD, & 
               maxSCF, thresh, max_diis, nBas, nOrb, nC, nO, nV, nR, Hc, ERI_AO, ERI_MO, ENuc, ERHF, eHF, cHF)

! Coupled-cluster module

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

  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)

! Local variables

  double precision              :: start_CC     ,end_CC       ,t_CC

!------------------------------------------------------------------------
! Perform CCD calculation
!------------------------------------------------------------------------

  if(doCCD) then

    call wall_time(start_CC)
    call CCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform DCD calculation
!------------------------------------------------------------------------

  if(doDCD) then

    call wall_time(start_CC)
    call DCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR, &
             ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for DCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform CCSD or CCSD(T) calculation
!------------------------------------------------------------------------

  if(doCCSDT) doCCSD = .true.

  if(doCCSD) then

    call wall_time(start_CC)
    call CCSD(dotest,maxSCF,thresh,max_diis,doCCSDT,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for CCSD or CCSD(T)= ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform direct ring CCD calculation
!------------------------------------------------------------------------

  if(dodrCCD) then

    call wall_time(start_CC)
    call drCCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for direct ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ring CCD calculation
!------------------------------------------------------------------------

  if(dorCCD) then

    call wall_time(start_CC)
    call rCCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for rCCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform crossed-ring CCD calculation
!------------------------------------------------------------------------

  if(docrCCD) then

    call wall_time(start_CC)
    call crCCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for crossed-ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ladder CCD calculation
!------------------------------------------------------------------------

  if(dolCCD) then

    call wall_time(start_CC)
    call lCCD(dotest,maxSCF,thresh,max_diis,nOrb,nC,nO,nV,nR,ERI_MO,ENuc,ERHF,eHF)
    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ladder CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform pair CCD calculation
!------------------------------------------------------------------------

  if(dopCCD) then

    call wall_time(start_CC)
    call pCCD(dotest, maxSCF, thresh, max_diis, nBas, nOrb, &
              nC, nO, nV, nR, Hc, ERI_AO, ENuc, ERHF, eHF, cHF)

    call wall_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for pair CCD = ',t_CC,' seconds'
    write(*,*)

  end if

end subroutine
