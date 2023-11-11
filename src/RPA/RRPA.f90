subroutine RRPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet, &  
                nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF,cHF,S)

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
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(dophRPA) then

    call wall_time(start_RPA)
    call phRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)
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
    call phRRPAx(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)
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
    call crRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ppRPA excitations
!------------------------------------------------------------------------

  if(doppRPA) then

    call wall_time(start_RPA)
    call ppRRPA(dotest,TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,EHF,ERI,dipole_int,epsHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

end subroutine
