subroutine complex_cRG0F2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,linearize,eta,regularize, &
                 nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,CAP,dipole_int,eHF)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

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

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  complex*16,intent(in)         :: eHF(nOrb)
  complex*16,intent(in)         :: CAP(nOrb,nOrb)
  complex*16,intent(in)         :: ERI(nOrb,nOrb,nOrb,nOrb)
  complex*16,intent(in)         :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  integer                       :: p
  double precision              :: Ec
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: Re_SigC(:)
  double precision,allocatable  :: Im_SigC(:)
  double precision,allocatable  :: Re_Z(:)
  double precision,allocatable  :: Im_Z(:)
  double precision,allocatable  :: Re_eGFlin(:)
  double precision, allocatable :: Im_eGFlin(:)
  double precision,allocatable  :: Re_eGF(:)
  double precision,allocatable  :: Im_eGF(:)
  double precision, allocatable :: Re_eHF(:)
  double precision, allocatable :: Im_eHF(:)

  ! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted G0F2 Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Memory allocation

  allocate(Re_SigC(nOrb),Im_SigC(nOrb), Re_Z(nOrb),Im_Z(nOrb),&
          Re_eGFlin(nOrb),Im_eGFlin(nOrb), Re_eGF(nOrb),Im_eGF(nOrb),Re_eHF(nOrb),Im_eHF(nOrb))
  Re_eHF(:) = real(eHF(:))
  Im_eHF(:) = aimag(eHF(:))

! Frequency-dependent second-order contribution

  call cRGF2_self_energy_diag(eta,nOrb,nC,nO,nV,nR,real(eHF),aimag(eHF),ERI,Re_SigC,Im_SigC,Re_Z,Im_Z)
  
  Re_eGFlin(:) = Re_eHF(:) + Re_Z(:)*Re_SigC(:) - Im_Z(:)*Im_SigC(:)
  Im_eGFlin(:) = Im_eHF(:) + Re_Z(:)*Im_SigC(:) + Im_Z(:)*Re_SigC(:)

  if(linearize) then

    write(*,*) '*** Quasiparticle energies obtained by linearization ***'

    Re_eGF(:) = Re_eGFlin(:)
    Im_eGF(:) = Im_eGFlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)
    write(*,*) "ONLY LINEARISATION IMPLEMENTED YET" 
    !call cRGF2_QP_graph(eta,nOrb,nC,nO,nV,nR,eHF,e_cap,ERI,Re_eGFlin,Im_eGFlin,eHF,e_cap,Re_eGF,Im_eGF,Re_Z,Im_Z)

  end if

  ! Print results

  call print_complex_cRG0F2(nOrb,nO,Re_eHF,Im_eHF,Re_SigC,Im_SigC,Re_eGF,Im_eGF,Re_Z,Im_Z,ENuc,ERHF,Ec)

  deallocate(Re_SigC,Im_SigC, Re_Z,Im_Z,&
          Re_eGFlin,Im_eGFlin, Re_eGF,Im_eGF)
end subroutine 
