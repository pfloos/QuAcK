subroutine complex_RRPA(use_gpu,dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,singlet,triplet,CVS,&
                nBas,nC,nO,nV,nR,nS,nCVS,ENuc,ERHF,ERI,dipole_int,CAP_MO,eHF,occupations)

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
  integer,intent(in)            :: nCVS
  integer,intent(in)            :: occupations(nO)
  double precision,intent(in)   :: ENuc
  complex*16,intent(in)         :: ERHF
  complex*16,intent(in)         :: eHF(nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: dipole_int(nBas,nBas,ncart)
  complex*16,intent(in)         :: CAP_MO(nBas,nBas)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

  if(CVS) then
    print *, "No CVS implemented"
  endif

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------

  if(dophRPA) then

    call wall_time(start_RPA)
    call complex_phRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,CAP_MO,eHF)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

end subroutine
