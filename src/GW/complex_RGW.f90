subroutine complex_RGW(dotest,docG0W0,doevGW,doqsGW,maxSCF,thresh,max_diis,doACFDT,                &
               exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_W,TDA,dBSE,dTDA,singlet,triplet, &
               linearize,eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,              &
               S,X,T,V,Hc,ERI_AO,ERI_MO,CAP_AO,CAP_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Restricted GW module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: docG0W0,doevGW,doqsGW

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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  complex*16,intent(in)         :: ERHF
  complex*16,intent(in)         :: eHF(nOrb)
  complex*16,intent(in)         :: cHF(nBas,nOrb)
  complex*16,intent(in)         :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: CAP_AO(nOrb,nOrb)
  complex*16,intent(in)         :: CAP_MO(nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  complex*16,intent(in)         :: dipole_int_MO(nOrb,nOrb,ncart)

! Local variables

  double precision              :: start_GW     ,end_GW       ,t_GW

  logical                       :: doccG0W0,doccGW

!------------------------------------------------------------------------
! Perform cG0W0 calculation
!------------------------------------------------------------------------

  if(docG0W0) then 
    call wall_time(start_GW)
    call complex_cRG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
               linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,CAP_MO,dipole_int_MO,eHF)
    call wall_time(end_GW)
  
    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for complex G0W0 = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGW calculation
!------------------------------------------------------------------------

  if(doevGW) then

    call wall_time(start_GW)
    call complex_evRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, &
               singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GW)

    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for evGW = ',t_GW,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGW calculation
!------------------------------------------------------------------------

  if(doqsGW) then 

    call wall_time(start_GW)
    call complex_qsRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2, &
               TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet,eta,doSRG,nNuc,ZNuc,rNuc,         &
               ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,ERI_MO,                  &
               dipole_int_AO,dipole_int_MO,PHF,cHF,eHF,                                      &
               CAP_AO,CAP_MO)
    call wall_time(end_GW)

    t_GW = end_GW - start_GW
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGW = ',t_GW,' seconds'
    write(*,*)

  end if

end subroutine
