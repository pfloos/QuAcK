subroutine GF(doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,unrestricted,renorm,maxSCF,thresh,max_diis,             &
              dophBSE,TDA,dBSE,dTDA,singlet,triplet,spin_conserved,spin_flip,linearize,eta,regularize,      & 
              nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
              dipole_int_AO,dipole_int,dipole_int_aa,dipole_int_bb,PHF,cHF,epsHF)

! Green's function module

  implicit none
  include 'parameters.h'

! Input variables

  logical                       :: doG0F2
  logical                       :: doevGF2
  logical                       :: doqsGF2
  logical                       :: doG0F3
  logical                       :: doevGF3
  logical                       :: unrestricted

  integer                       :: renorm
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)

  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  double precision              :: start_GF     ,end_GF       ,t_GF

!------------------------------------------------------------------------
! Compute G0F2 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F2) then

    call cpu_time(start_GF)
    if(unrestricted) then
      call UG0F2(dophBSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,linearize,eta,regularize, &
                 nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_aaaa,ERI_aabb,ERI_bbbb,  & 
                 dipole_int_aa,dipole_int_bb,epsHF)
    else
      call G0F2(dophBSE,TDA,dBSE,dTDA,singlet,triplet,linearize,eta,regularize, & 
                nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,epsHF)
    end if
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF2 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF2) then

    call cpu_time(start_GF)
    if(unrestricted) then
      call evUGF2(maxSCF,thresh,max_diis,dophBSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,  &
                  eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
                  dipole_int_aa,dipole_int_bb,cHF,epsHF)
    else
      call evGF2(dophBSE,TDA,dBSE,dTDA,maxSCF,thresh,max_diis, & 
                 singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF, & 
                 ERI,dipole_int,epsHF)
    end if
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGF2 calculation
!------------------------------------------------------------------------

  if(doqsGF2) then 

    call cpu_time(start_GF)
    if(unrestricted) then
      call qsUGF2(maxSCF,thresh,max_diis,dophBSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,regularize, &
                  nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,                        & 
                  ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,epsHF)
    else
      call qsGF2(maxSCF,thresh,max_diis,dophBSE,TDA,dBSE,dTDA,singlet,triplet,eta,regularize,nNuc,ZNuc,rNuc,ENuc, & 
                 nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI,dipole_int_AO,dipole_int,PHF,cHF,epsHF)
    end if
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute G0F3 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F3) then

    call cpu_time(start_GF)
    if(unrestricted) then
      print*,'Unrestricted version of G0F3 not yet implemented! Sorry.'
    else
      call G0F3(renorm,nBas,nC,nO,nV,nR,ERI,epsHF)
    end if
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF3 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF3 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF3) then

    call cpu_time(start_GF)
    if(unrestricted) then
      print*,'Unrestricted version of evGF3 not yet implemented! Sorry.'
    else
      call evGF3(maxSCF,thresh,max_diis,renorm,nBas,nC,nO,nV,nR,ERI,epsHF)
    end if
    call cpu_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GF3 = ',t_GF,' seconds'
    write(*,*)

  end if

end subroutine
