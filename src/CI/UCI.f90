subroutine UCI(dotest,doCIS,doCIS_D,doCID,doCISD,doFCI,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS, & 
               ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,epsHF,EHF,cHF,S,F)

! Configuration interaction module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doCIS
  logical,intent(in)            :: doCIS_D
  logical,intent(in)            :: doCID
  logical,intent(in)            :: doCISD
  logical,intent(in)            :: doFCI

  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: F(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart,nspin)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart,nspin)

! Local variables

  double precision              :: start_CI     ,end_CI       ,t_CI

!------------------------------------------------------------------------
! Compute CIS excitations
!------------------------------------------------------------------------

  if(doCIS) then

    call wall_time(start_CI)
    call UCIS(dotest,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb, & 
              ERI_bbbb,dipole_int_aa,dipole_int_bb,epsHF,cHF,S)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CID excitations
!------------------------------------------------------------------------

  if(doCID) then

    call wall_time(start_CI)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CID = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute CISD excitations
!------------------------------------------------------------------------

  if(doCISD) then

    call wall_time(start_CI)
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CISD = ',t_CI,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute FCI 
!------------------------------------------------------------------------

  if(doFCI) then

    call wall_time(start_CI)
    write(*,*) ' FCI is not yet implemented! Sorry.'
    call wall_time(end_CI)

    t_CI = end_CI - start_CI
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for FCI = ',t_CI,' seconds'
    write(*,*)

  end if

end subroutine
