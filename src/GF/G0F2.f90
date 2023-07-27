subroutine G0F2(dophBSE,doppBSE,TDA,dBSE,dTDA,singlet,triplet,linearize,eta,regularize, &
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|     One-shot second-order Green function     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(SigC(nBas),Z(nBas),eGF(nBas))

  if(linearize) then 
  
     write(*,*) '*** Quasiparticle equation will be linearized ***'
     write(*,*)

  end  if

! Frequency-dependent second-order contribution

  if(regularize) then 

    call regularized_self_energy_GF2_diag(eta,nBas,nC,nO,nV,nR,eHF,eHF,ERI,SigC,Z)

  else

    call GF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,eHF,eHF,ERI,SigC,Z)

  end if
  
  if(linearize) then

    eGF(:) = eHF(:) + Z(:)*SigC(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
    write(*,*)

    call GF2_QP_graph(eta,nBas,nC,nO,nV,nR,eHF,ERI,eGF)

  end if

  ! Print results

  call MP2(regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,eGF,Ec)
  call print_G0F2(nBas,nO,eHF,SigC,eGF,Z,ENuc,ERHF,Ec)

! Perform BSE2 calculation

  if(dophBSE) then 
  
    call GF2_phBSE2(TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  correlation energy           =',sum(EcBSE(:))
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@G0F2  total energy                 =',ENuc + EHF + sum(EcBSE(:))
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Perform ppBSE2 calculation

  if(doppBSE) then 
   
    call GF2_ppBSE2(TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy (triplet) =',3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 correlation energy           =',EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0F2 total energy                 =',ENuc + ERHF + EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine 
