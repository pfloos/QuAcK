subroutine UGW(dotest,doG0W0,doevGW,doqsGW,doufG0W0,doufGW,maxSCF,thresh,max_diis,doACFDT,                 &
              exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
              linearize,eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,                 & 
              ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)

! GW module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doG0W0
  logical,intent(in)            :: doevGW
  logical,intent(in)            :: doqsGW
  logical,intent(in)            :: doufG0W0
  logical,intent(in)            :: doufGW

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

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

  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  double precision              :: start_GW     ,end_GW       ,t_GW

!------------------------------------------------------------------------
! Perform G0W0 calculatiom
!------------------------------------------------------------------------

  if(doG0W0) then
    
    call wall_time(start_GW)
    call UG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip, & 
               linearize,eta,doSRG,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,S,ERI_aaaa,ERI_aabb,ERI_bbbb,            & 
               dipole_int_aa,dipole_int_bb,cHF,eHF)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0W0 = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGW calculation
!------------------------------------------------------------------------

  if(doevGW) then

    call wall_time(start_GW)
    call evUGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_W,TDA,dBSE,dTDA, & 
              spin_conserved,spin_flip,linearize,eta,doSRG,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,S,  & 
              ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eHF)
    call wall_time(end_GW)

    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGW = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGW calculation
!------------------------------------------------------------------------

  if(doqsGW) then 

    call wall_time(start_GW)
    call qsUGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_W,TDA,dBSE,dTDA, & 
               spin_conserved,spin_flip,eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,         & 
               S,X,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)
    call wall_time(end_GW)

    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGW = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ufG0W0 calculatiom
!------------------------------------------------------------------------

  if(doufG0W0) then
    
    call wall_time(start_GW)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ufG0W0 = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ufGW calculatiom
!------------------------------------------------------------------------

  if(doufGW) then
    
    call wall_time(start_GW)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for ufGW = ',t_GW,' seconds'
    write(*,*)

  end if

end subroutine
