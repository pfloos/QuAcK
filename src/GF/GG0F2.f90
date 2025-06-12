subroutine GG0F2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,linearize,eta,regularize, &
                 nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

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
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: Ec
  double precision              :: EcBSE
  double precision,allocatable  :: eGFlin(:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)

! Hello world


  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Generalized G0F2 Calculation *'
  write(*,*)'********************************'
  write(*,*)

! Memory allocation

  allocate(SigC(nBas),Z(nBas),eGFlin(nBas),eGF(nBas))

! Frequency-dependent second-order contribution

  if(regularize) then 

!   call GF2_reg_self_energy_diag(eta,nBas,nC,nO,nV,nR,eHF,ERI,SigC,Z)

  else

    call GGF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,eHF,ERI,SigC,Z)

  end if
  
  eGFlin(:) = eHF(:) + Z(:)*SigC(:)

  if(linearize) then

    write(*,*) '*** Quasiparticle energies obtained by linearization ***'

    eGF(:) = eGFlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)

    call GGF2_QP_graph(eta,nBas,nC,nO,nV,nR,eHF,ERI,eGFlin,eHF,eGF,Z)

  end if

  ! Print results

  call GMP2(.false.,regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,eGF,Ec)
  call print_GG0F2(nBas,nC,nO,nV,nR,eHF,SigC,eGF,Z,ENuc,EGHF,Ec)

! Perform BSE@GF2 calculation

  if(dophBSE) then 
  
    call GGF2_phBSE(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@GG0F2  correlation energy          =',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@GG0F2  total energy                =',ENuc + EHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Perform ppBSE@GF2 calculation

  if(doppBSE) then 
   
    call GGF2_ppBSE(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@GG0F2 correlation energy          =',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@GG0F2 total energy                =',ENuc + EGHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('G','G0F2 correlation energy',Ec)
    call dump_test_value('G','G0F2 HOMO energy',eGF(nO))
    call dump_test_value('G','G0F2 LUMO energy',eGF(nO+1))

  end if

end subroutine 
