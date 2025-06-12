subroutine GGT(dotest,doG0T0pp,doevGTpp,doqsGTpp,doG0T0eh,maxSCF,thresh,max_diis, & 
               doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,            &  
               TDA_T,TDA,dBSE,dTDA,linearize,eta,doSRG,nNuc,ZNuc,rNuc,ENuc,       & 
               nBas,nBas2,nC,nO,nV,nR,nS,EGHF,S,X,T,V,Hc,                         & 
               ERI_AO,ERI,dipole_int_AO,dipole_int,PHF,cHF,eHF)
  
! GT module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doG0T0pp
  logical,intent(in)            :: doevGTpp
  logical,intent(in)            :: doqsGTpp

  logical,intent(in)            :: doG0T0eh

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas2)
  double precision,intent(in)   :: cHF(nBas2,nBas2)
  double precision,intent(in)   :: PHF(nBas2,nBas2)
  double precision,intent(in)   :: S(nBas2,nBas2)
  double precision,intent(in)   :: T(nBas2,nBas2)
  double precision,intent(in)   :: V(nBas2,nBas2)
  double precision,intent(in)   :: Hc(nBas2,nBas2)
  double precision,intent(in)   :: X(nBas2,nBas2)
  double precision,intent(in)   :: ERI_AO(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: ERI(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: dipole_int_AO(nBas2,nBas2,ncart)
  double precision,intent(in)   :: dipole_int(nBas2,nBas2,ncart)

! Local variables

  double precision              :: start_GT     ,end_GT       ,t_GT

!------------------------------------------------------------------------
! Perform G0T0pp calculatiom
!------------------------------------------------------------------------

  if(doG0T0pp) then
    call wall_time(start_GT)
    call GG0T0pp(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, &
              linearize,eta,doSRG,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0T0pp = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGTpp calculation
!------------------------------------------------------------------------

  ! if(doevGTpp) then

  !   call wall_time(start_GT)
  !   call evGGTpp(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, &
  !              linearize,eta,doSRG,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
  !   call wall_time(end_GT)

  !   t_GT = end_GT - start_GT
  !   write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGT = ',t_GT,' seconds'
  !   write(*,*)

  ! end if

!------------------------------------------------------------------------
! Perform qsGTpp calculation
!------------------------------------------------------------------------

  ! if(doqsGTpp) then 

  !   call wall_time(start_GT)
  !   call qsGGTpp(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, & 
  !              eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nBas2,nC,nO,nV,nR,nS,EGHF,S,X,T,V,Hc,ERI_AO,ERI,                       & 
  !              dipole_int_AO,dipole_int,PHF,cHF,eHF)
  !   call wall_time(end_GT)

  !   t_GT = end_GT - start_GT
  !   write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGT = ',t_GT,' seconds'
  !   write(*,*)

  ! end if


!------------------------------------------------------------------------
! Perform G0T0eh calculatiom
!------------------------------------------------------------------------

  if(doG0T0eh) then
    call wall_time(start_GT)
    call GG0T0eh(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, &
              linearize,eta,doSRG,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0T0eh = ',t_GT,' seconds'
    write(*,*)

  end if

end subroutine
