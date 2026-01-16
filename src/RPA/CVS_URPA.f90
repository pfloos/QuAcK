subroutine CVS_URPA(dotest,dophRPA,dophRPAx,docrRPA,doppRPA,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip, &  
                nBas,nC,nO,nV,nR,nS,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,cHF,S,occupations)

! MERGE THIS WITH THE URPA BRANCH

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
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  double precision              :: start_RPA    ,end_RPA      ,t_RPA

!------------------------------------------------------------------------
! Compute (direct) RPA excitations
!------------------------------------------------------------------------
print *, "Hello world"
  if(dophRPA) then

    call wall_time(start_RPA)
    call CVS_phURPA(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,ENuc,EUHF, &
                ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,cHF,S,occupations)
    call wall_time(end_RPA)

    t_RPA = end_RPA - start_RPA
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPA = ',t_RPA,' seconds'
    write(*,*)

  end if

!1!------------------------------------------------------------------------
!1! Compute RPAx (RPA with exchange) excitations
!1!------------------------------------------------------------------------
!1
!1  if(dophRPAx) then
!1
!1    call wall_time(start_RPA)
!1    call phURPAx(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,ENuc,EUHF, &
!1                 ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,cHF,S)
!1    call wall_time(end_RPA)
!1
!1    t_RPA = end_RPA - start_RPA
!1    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RPAx = ',t_RPA,' seconds'
!1    write(*,*)
!1
!1  end if
!1
!1!------------------------------------------------------------------------
!1! Compute crRPA excitations
!1!------------------------------------------------------------------------
!1
!1  if(docrRPA) then
!1
!1    call wall_time(start_RPA)
!1    write(*,*) 'Unrestricted version of crRPA not yet implemented! Sorry.'
!1    call wall_time(end_RPA)
!1
!1    t_RPA = end_RPA - start_RPA
!1    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
!1    write(*,*)
!1
!1  end if
!1
!1!------------------------------------------------------------------------
!1! Compute ppRPA excitations
!1!------------------------------------------------------------------------
!1
!1  if(doppRPA) then
!1
!1    call wall_time(start_RPA)
!1    call ppURPA(dotest,TDA,doACFDT,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF)
!1    call wall_time(end_RPA)
!1
!1    t_RPA = end_RPA - start_RPA
!1    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for pp-RPA = ',t_RPA,' seconds'
!1    write(*,*)
!1
!1  end if

end subroutine
