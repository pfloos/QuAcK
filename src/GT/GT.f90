subroutine GT(doG0T0pp,doevGTpp,doqsGTpp,doG0T0eh,doevGTeh,doqsGTeh,unrestricted,maxSCF,thresh,max_diis,doACFDT,  &
              exchange_kernel,doXBS,dophBSE,dophBSE2,doppBSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet,spin_conserved,spin_flip, &
              linearize,eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,               &
              ERI_AO,ERI,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int,dipole_int_aa,dipole_int_bb,    &
              PHF,cHF,epsHF)

! T-matrix module

  implicit none
  include 'parameters.h'

! Input variables

  logical                       :: doG0T0pp
  logical                       :: doevGTpp
  logical                       :: doqsGTpp
  logical                       :: doG0T0eh
  logical                       :: doevGTeh
  logical                       :: doqsGTeh
  logical                       :: unrestricted

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
  double precision,intent(inout):: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)


! Local variables 

  double precision              :: start_GT     ,end_GT       ,t_GT

!------------------------------------------------------------------------
! Perform G0T0pp calculatiom
!------------------------------------------------------------------------

  if(doG0T0pp) then
    
    call wall_time(start_GT)
    if(unrestricted) then 
       call UG0T0pp(doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip,      &
                    linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb, &
                    dipole_int_aa,dipole_int_bb,PHF,cHF,epsHF)
    else
      call G0T0pp(doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
                  linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI,dipole_int,PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0T0 = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGTpp calculatiom
!------------------------------------------------------------------------

  if(doevGTpp) then
    
    call wall_time(start_GT)
    if(unrestricted) then
      call evUGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
                   eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb, & 
                   PHF,cHF,epsHF)
    else
      call evGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet, &
                  eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI,dipole_int,PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGTpp calculation
!------------------------------------------------------------------------

  if(doqsGTpp) then 

    call wall_time(start_GT)
    if(unrestricted) then
      call qsUGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip, & 
                   eta,regularize,nBas,nC,nO,nV,nR,nS,nNuc,ZNuc,rNuc,ENuc,EHF,S,X,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,   &
                   dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,epsHF) 
    else
      call qsGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet,          &
                  eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI,dipole_int_AO,dipole_int, & 
                  PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for qsGT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform G0T0eh calculatiom
!------------------------------------------------------------------------

  if(doG0T0eh) then
    
    call wall_time(start_GT)
    if(unrestricted) then 
      print*,'Unrestricted version of G0T0eh not yet implemented! Sorry.'
    else
      call G0T0eh(doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, &
                  linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI,dipole_int,PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)
  
    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for G0T0 = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform evGTeh calculation
!------------------------------------------------------------------------

  if(doevGTeh) then

    call wall_time(start_GT)
    if(unrestricted) then 
      print*,'Unrestricted version of evGTeh not yet implemented! Sorry.'
    else
      call evGTeh(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, &
                  singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI_AO,ERI,dipole_int,       & 
                  PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for evGT = ',t_GT,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGTeh calculation
!------------------------------------------------------------------------

  if(doqsGTeh) then 

    call wall_time(start_GT)
    if(unrestricted) then 
      print*,'Unrestricted version of qsGTeh not yet implemented! Sorry.'
    else
      call qsGTeh(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,singlet,triplet, &
                  eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EHF,S,X,T,V,Hc,ERI_AO,ERI,dipole_int_AO,dipole_int, & 
                  PHF,cHF,epsHF)
    end if
    call wall_time(end_GT)

    t_GT = end_GT - start_GT
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGW = ',t_GT,' seconds'
    write(*,*)

  end if

end subroutine
