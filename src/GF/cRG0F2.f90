subroutine cRG0F2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,linearize,eta,regularize, &
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
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)
  double precision,intent(in)   :: CAP(nOrb,nOrb)

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
  double precision, allocatable :: e_CAP(:)

  ! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted G0F2 Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Memory allocation

  allocate(Re_SigC(nOrb),Im_SigC(nOrb), Re_Z(nOrb),Im_Z(nOrb),&
          Re_eGFlin(nOrb),Im_eGFlin(nOrb), Re_eGF(nOrb),Im_eGF(nOrb),e_CAP(nOrb))
  do p = 1, nOrb
        e_CAP(p) = CAP(p,p)
        write(*,*) p, e_CAP(p)
  end do
  
! Frequency-dependent second-order contribution

  if(regularize) then 
   write(*,*) "Regularisation not implemented (yet)"
   ! call RGF2_reg_self_energy_diag(eta,nOrb,nC,nO,nV,nR,eHF,ERI,SigC,Z)

  else

    call cRGF2_self_energy_diag(eta,nOrb,nC,nO,nV,nR,eHF,ERI,Re_SigC,Im_SigC,Re_Z,Im_Z,e_CAP)

  end if
  
  Re_eGFlin(:) = eHF(:) + Re_Z(:)*Re_SigC(:) - Im_Z(:)*Im_SigC(:)
  Im_eGFlin(:) = e_CAP(:) + Re_Z(:)*Im_SigC(:) + Im_Z(:)*Re_SigC(:)

  if(linearize) then

    write(*,*) '*** Quasiparticle energies obtained by linearization ***'

    Re_eGF(:) = Re_eGFlin(:)
    Im_eGF(:) = Im_eGFlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** NOT IMPLEMEMTED YET'
    write(*,*)
    
    !call RGF2_QP_graph(eta,nOrb,nC,nO,nV,nR,eHF,ERI,eGFlin,eHF,eGF,Z)

  end if

  ! Print results

!  call RMP2(.false.,regularize,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eGF,Ec)
  call print_cRG0F2(nOrb,nO,eHF,e_CAP,Re_SigC,Im_SigC,Re_eGF,Im_eGF,Re_Z,Im_Z,ENuc,ERHF,Ec)

! Perform BSE@GF2 calculation
!
!  if(dophBSE) then 
!  
!    call RGF2_phBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy (singlet) =',EcBSE(1)
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy (triplet) =',EcBSE(2)
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy           =',sum(EcBSE)
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  total energy                 =',ENuc + ERHF + sum(EcBSE)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Perform ppBSE@GF2 calculation
!
!  if(doppBSE) then 
!   
!     call RGF2_ppBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)
!
!    EcBSE(2) = 3d0*EcBSE(2)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy (singlet) =',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy (triplet) =',EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy           =',sum(EcBSE),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 total energy                 =',ENuc + ERHF + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Testing zone
!
!  if(dotest) then
!
!    call dump_test_value('R','G0F2 correlation energy',Ec)
!    call dump_test_value('R','G0F2 HOMO energy',eGF(nO))
!    call dump_test_value('R','G0F2 LUMO energy',eGF(nO+1))
!
!  end if
!
!  deallocate(SigC, Z, eGFlin, eGF)
!
end subroutine 
