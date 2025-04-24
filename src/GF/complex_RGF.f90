subroutine complex_RGF(dotest,docG0F2,doevGF2,doqsGF2,maxSCF,                           &
               thresh,max_diis,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,linearize, &
               eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,        &
               S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF,        &
               CAP_AO,CAP_MO)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest
  logical,intent(in)            :: docG0F2,doevGF2,doqsGF2

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis

  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

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
  complex*16,intent(in)         :: S(nBas,nBas)
  complex*16,intent(in)         :: CAP_AO(nBas,nBas)
  complex*16,intent(in)         :: CAP_MO(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  complex*16,intent(in)         :: dipole_int_MO(nOrb,nOrb,ncart)

! Local variables

  double precision              :: start_GF     ,end_GF       ,t_GF


!------------------------------------------------------------------------
! Compute complex G0F2 electronic binding energies
!------------------------------------------------------------------------

  if(docG0F2) then

    call wall_time(start_GF)
    call complex_cRG0F2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet, &
               linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,    &
               ENuc,ERHF,ERI_MO,CAP_MO,dipole_int_MO,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if
  if(doevGF2) then

    call wall_time(start_GF)
    call complex_evRGF2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,maxSCF,thresh,max_diis,singlet,triplet, &
                 linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

  if(doqsGF2) then

    call wall_time(start_GF)
    call complex_qsRGF2(dotest,maxSCF,thresh,max_diis,dophBSE,doppBSE,TDA,  &
                  dBSE,dTDA,singlet,triplet,eta,regularize,nNuc,ZNuc, &
                  rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc, & 
                  ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF, &
                  CAP_AO,CAP_MO)
    call wall_time(end_GF)

    t_GF = end_GF - start_GF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for GF2 = ',t_GF,' seconds'
    write(*,*)

  end if

end subroutine
