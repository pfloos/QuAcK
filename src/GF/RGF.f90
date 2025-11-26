subroutine RGF(dotest,doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3,renorm,maxSCF,              &
               thresh,max_diis,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,linearize, &
               eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,             &
               S,X,T,V,Hc,ERI_AO,ERI_MO,CAP,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Green's function module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doG0F2
  logical,intent(in)            :: doevGF2
  logical,intent(in)            :: doqsGF2
  logical,intent(in)            :: doG0F3
  logical,intent(in)            :: doevGF3

  integer                       :: renorm
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
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

  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: CAP(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_MO(nOrb,nOrb,ncart)

! Local variables

  double precision              :: start_GF     ,end_GF       ,t_GF

!------------------------------------------------------------------------
! Compute G0F2 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F2) then

    call wall_time(start_GF)
    call RG0F2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet, &
               linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,         &
               ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF2 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF2) then

    call wall_time(start_GF)
    call evRGF2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,maxSCF,thresh,max_diis,  & 
                singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS, &
                ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Perform qsGF2 calculation
!------------------------------------------------------------------------

  if(doqsGF2) then 

    call wall_time(start_GF)
    call qsRGF2(dotest,maxSCF,thresh,max_diis,dophBSE,doppBSE,TDA,  &
                dBSE,dTDA,singlet,triplet,eta,doSRG,nNuc,ZNuc,      &
                rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc, & 
                ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for qsGF2 = ',t_GF,' seconds'
    write(*,*)

 end if

!------------------------------------------------------------------------
! Compute G0F3 electronic binding energies
!------------------------------------------------------------------------

  if(doG0F3) then

    call wall_time(start_GF)
    call RG0F3(dotest,renorm,nBas,nOrb,nC,nO,nV,nR,ERI_MO, eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF3 = ',t_GF,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute evGF3 electronic binding energies
!------------------------------------------------------------------------

  if(doevGF3) then

    call wall_time(start_GF)
    call evRGF3(dotest,maxSCF,thresh,max_diis,renorm,nBas,nOrb,nC,nO,nV,nR,ERI_MO,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF3 = ',t_GF,' seconds'
    write(*,*)

  end if

end subroutine
