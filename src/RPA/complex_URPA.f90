subroutine complex_URPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,CVS, &  
                nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,&
                dipole_int_aa,dipole_int_bb,CAP_MO,eHF,cHF,S,occupations)

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
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: CVS
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  complex*16,intent(in)         :: EUHF
  complex*16,intent(in)         :: eHF(nBas,nspin)
  complex*16,intent(in)         :: cHF(nBas,nBas,nspin)
  complex*16,intent(in)         :: CAP_MO(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  complex*16 ,intent(in)        :: ERI_aaaa(nBas,nBas,nBas,nBas)
  complex*16 ,intent(in)        :: ERI_aabb(nBas,nBas,nBas,nBas)
  complex*16 ,intent(in)        :: ERI_bbbb(nBas,nBas,nBas,nBas)
  complex*16 ,intent(in)        :: dipole_int_aa(nBas,nBas,ncart)
  complex*16 ,intent(in)        :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------
  
  if(dophRPA) then

    call wall_time(start_RPA)
    call complex_phURPA(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF, &
                ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,CAP_MO,eHF,cHF,S,occupations)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute RPAx (RPA with exchange) excitations
!------------------------------------------------------------------------

  if(dophRPAx .and. CVS) then
    print *, 'No phRPAx yet for CAP/complex'
   ! call wall_time(start_RPA)
   ! call CVS_phURPAx(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF, &
   !             ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,cHF,S,occupations)
   ! call wall_time(end_RPA)

   ! t_RPA = end_RPA - start_RPA
   ! write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPAx = ',t_RPA,' seconds'
   ! write(*,*)

  end if

end subroutine
