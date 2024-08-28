
! ---

subroutine RGT(dotest, doG0T0pp, doevGTpp, doqsGTpp, doufG0T0pp, doG0T0eh, doevGTeh, doqsGTeh, &
               maxSCF, thresh, max_diis, doACFDT, exchange_kernel, doXBS, dophBSE, dophBSE2,   &
               doppBSE, TDA_T, TDA, dBSE, dTDA, singlet, triplet, linearize, eta, regularize,  &
               nNuc, ZNuc, rNuc, ENuc, nBas_AOs, nBas_MOs, nC, nO, nV, nR, nS, ERHF, S, X, T,  &
               V, Hc, ERI_AO, ERI_MO, dipole_int_AO, dipole_int_MO, PHF, cHF, eHF)

! T-matrix module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doG0T0pp
  logical,intent(in)            :: doevGTpp
  logical,intent(in)            :: doqsGTpp
  logical,intent(in)            :: doufG0T0pp
  logical,intent(in)            :: doG0T0eh
  logical,intent(in)            :: doevGTeh
  logical,intent(in)            :: doqsGTeh

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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas_MOs)
  double precision,intent(in)   :: cHF(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: PHF(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: T(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: V(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: ERI_AO(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: ERI_MO(nBas_MOs,nBas_MOs,nBas_MOs,nBas_MOs)
  double precision,intent(in)   :: dipole_int_AO(nBas_AOs,nBas_AOs,ncart)
  double precision,intent(in)   :: dipole_int_MO(nBas_MOs,nBas_MOs,ncart)


! Local variables 

  double precision              :: start_GT     ,end_GT       ,t_GT

!------------------------------------------------------------------------
! Perform G0T0pp calculatiom
!------------------------------------------------------------------------

  if(doG0T0pp) then
    
    call wall_time(start_GT)
    call RG0T0pp(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
                linearize,eta,regularize,nBas_MOs,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for G0T0pp = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGTpp calculatiom
!------------------------------------------------------------------------

  if(doevGTpp) then
    
    call wall_time(start_GT)
    call evRGTpp(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet, &
                 linearize,eta,regularize,nBas_MOs,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for evGTpp = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGTpp calculation
!------------------------------------------------------------------------

  if(doqsGTpp) then 

    call wall_time(start_GT)
    call qsRGTpp(dotest, maxSCF, thresh, max_diis, doACFDT, exchange_kernel, doXBS, dophBSE, TDA_T, TDA, dBSE, &
                 dTDA, singlet, triplet, eta, regularize, nNuc, ZNuc, rNuc, ENuc, nBas_AOs, nBas_MOs, nC, nO,  &
                 nV, nR, nS, ERHF, S, X, T, V, Hc, ERI_AO, ERI_MO, dipole_int_AO, dipole_int_MO, PHF, cHF, eHF)
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGTpp = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform ufG0T0pp calculatiom
!------------------------------------------------------------------------

  if(doufG0T0pp) then
    
    call wall_time(start_GT)
    call ufG0T0pp(dotest,TDA_T,nBas_MOs,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ufG0T0pp = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform G0T0eh calculatiom
!------------------------------------------------------------------------

  if(doG0T0eh) then
    
    call wall_time(start_GT)
    call RG0T0eh(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
                 linearize,eta,regularize,nBas_MOs,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for G0T0eh = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGTeh calculation
!------------------------------------------------------------------------

  if(doevGTeh) then

    call wall_time(start_GT)
    call evRGTeh(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, &
                 singlet,triplet,linearize,eta,regularize,nBas_MOs,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for evGTeh = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGTeh calculation
!------------------------------------------------------------------------

  if(doqsGTeh) then 

    call wall_time(start_GT)
    call qsRGTeh(dotest, maxSCF, thresh, max_diis, doACFDT, exchange_kernel, doXBS, dophBSE, &
                 dophBSE2, TDA_T, TDA, dBSE, dTDA, singlet, triplet, eta, regularize, nNuc,  &
                 ZNuc, rNuc, ENuc, nBas_AOs, nBas_MOs, nC, nO, nV, nR, nS, ERHF, S, X, T, V, &
                 Hc, ERI_AO, ERI_MO, dipole_int_AO, dipole_int_MO, PHF, cHF, eHF)
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGTeh = ',t_GT,' seconds'
    write(*,*)

  end if

end subroutine
