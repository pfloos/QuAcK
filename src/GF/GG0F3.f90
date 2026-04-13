subroutine GG0F3(dotest,linearize,eta,doSRG,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! Perform a one-shot third-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)

! Local variables

  double precision              :: Ec
  double precision              :: flow
  double precision,allocatable  :: eGFlin(:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Generalized G0F3 Calculation *'
  write(*,*)'********************************'
  write(*,*)

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G0F3 scheme ***'
    write(*,*)

  end if
  
! Memory allocation

  allocate(SigC(nBas), Z(nBas), eGFlin(nBas), eGF(nBas))
  SigC(:) = 0d0
  Z(:) = 0d0
  eGFlin(:) = 0d0
  eGF(:) = 0d0

! Frequency-dependent third-order contribution

  if(doSRG) then 

    ! TODO call GGF3_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,eHF,ERI,Ec,SigC,Z)

  else

    call GGF3_self_energy_diag(eta,nBas,nC,nO,nV,nR,eHF,ERI,SigC,Z)

  end if
  
  eGFlin(:) = eHF(:) + Z(:)*SigC(:)

  if(linearize) then

    write(*,*) '*** Quasiparticle energies obtained by linearization ***'

    eGF(:) = eGFlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)

    call GGF3_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,eHF,ERI,eGFlin,eHF,eGF,Z)

  end if

  ! Print results

  call print_GG0F3(nBas,nC,nO,nV,nR,eHF,SigC,eGF,Z,ENuc,EGHF,Ec)

! Testing zone

  if(dotest) then

    call dump_test_value('R','G0F2 correlation energy',Ec)
    call dump_test_value('R','G0F2 HOMO energy',eGF(nO))
    call dump_test_value('R','G0F2 LUMO energy',eGF(nO+1))

  end if

  deallocate(SigC, Z, eGFlin, eGF)
  
end subroutine GG0F3
