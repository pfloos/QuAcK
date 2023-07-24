subroutine CC(doCCD,dopCCD,doDCD,doCCSD,doCCSDT,do_drCCD,do_rCCD,do_crCCD,do_lCCD, & 
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
  logical                       :: do_drCCD
  logical                       :: do_rCCD
  logical                       :: do_crCCD
  logical                       :: do_lCCD

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
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

    call cpu_time(start_CC)
    call CCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform DCD calculation
!------------------------------------------------------------------------

  if(doDCD) then

    call cpu_time(start_CC)
    call DCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR, & 
             ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for DCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform CCSD or CCSD(T) calculation
!------------------------------------------------------------------------

  if(doCCSDT) doCCSD = .true.

  if(doCCSD) then

    call cpu_time(start_CC)
    call CCSD(maxSCF,thresh,max_diis,doCCSDT,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCSD or CCSD(T)= ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform direct ring CCD calculation
!------------------------------------------------------------------------

  if(do_drCCD) then

    call cpu_time(start_CC)
    call drCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for direct ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ring CCD calculation
!------------------------------------------------------------------------

  if(do_rCCD) then

    call cpu_time(start_CC)
    call rCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for rCCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform crossed-ring CCD calculation
!------------------------------------------------------------------------

  if(do_crCCD) then

    call cpu_time(start_CC)
    call crCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for crossed-ring CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ladder CCD calculation
!------------------------------------------------------------------------

  if(do_lCCD) then

    call cpu_time(start_CC)
    call lCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ladder CCD = ',t_CC,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform pair CCD calculation
!------------------------------------------------------------------------

  if(dopCCD) then

    call cpu_time(start_CC)
    call pCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call cpu_time(end_CC)

    t_CC = end_CC - start_CC
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pair CCD = ',t_CC,' seconds'
    write(*,*)

  end if

end subroutine