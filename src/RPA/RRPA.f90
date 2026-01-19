subroutine RRPA(use_gpu,dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet,CVS,&
                nBas,nC,nO,nV,nR,nS,nCVS,ENuc,ERHF,ERI,dipole_int,eHF,occupations)

! Random-phase approximation module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: use_gpu

  logical,intent(in)            :: dotest

  logical,intent(in)            :: dophRPA
  logical,intent(in)            :: dophRPAx
  logical,intent(in)            :: docrRPA
  logical,intent(in)            :: doppRPA

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: CVS
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: occupations(nO,nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

  if(CVS) then
    print *, "CVS for restricted branch is not implemented yet, sry... MOM and the RPA methods are not available."
    stop
  endif

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(dophRPA) then

    call wall_time(start_RPA)
    if (use_gpu) then
      call phRRPA_GPU(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)
    else
      call phRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)
    endif
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute RPAx (RPA with exchange) excitations
!------------------------------------------------------------------------

  if(dophRPAx) then

    call wall_time(start_RPA)
    call phRRPAx(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPAx = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute crRPA excitations
!------------------------------------------------------------------------

  if(docrRPA) then

    call wall_time(start_RPA)
    call crRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ppRPA excitations
!------------------------------------------------------------------------

  if(doppRPA) then

    call wall_time(start_RPA)
    call ppRRPA(dotest,TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,dipole_int,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

end subroutine
