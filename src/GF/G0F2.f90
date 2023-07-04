subroutine G0F2(BSE,TDA,dBSE,dTDA,evDyn,singlet,triplet,linearize,eta,regularize, &
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
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
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision              :: Ec
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: eGF2lin(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|     One-shot second-order Green function     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(SigC(nBas),Z(nBas),eGF2(nBas),eGF2lin(nBas))

  if(linearize) then 
  
     write(*,*) '*** Quasiparticle equation will be linearized ***'
     write(*,*)

  end  if

! Frequency-dependent second-order contribution

  if(regularize) then 

    call regularized_self_energy_GF2_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,eHF,ERI,SigC,Z)

  else

    call GF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,eHF,ERI,SigC,Z)

  end if

  eGF2lin(:) = eHF(:) + Z(:)*SigC(:)
  
  if(linearize) then

    eGF2(:) = eGF2lin(:)

 else

    write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
    write(*,*)

    call QP_graph_GF2(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2lin,ERI,eGF2)


  end if

  ! Print results

  call MP2(regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,eGF2,Ec)
  call print_G0F2(nBas,nO,eHF,SigC,eGF2,Z,ENuc,ERHF,Ec)

! Perform BSE2 calculation

  if(BSE) then

    call GF2_phBSE2(TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF2,EcBSE)

  end if

end subroutine 
