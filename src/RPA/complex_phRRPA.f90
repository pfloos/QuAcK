subroutine complex_phRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform a direct random phase approximation calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)           :: dotest

  logical,intent(in)           :: TDA
  logical,intent(in)           :: doACFDT
  logical,intent(in)           :: exchange_kernel
  logical,intent(in)           :: singlet
  logical,intent(in)           :: triplet
  integer,intent(in)           :: nBas
  integer,intent(in)           :: nC
  integer,intent(in)           :: nO
  integer,intent(in)           :: nV
  integer,intent(in)           :: nR
  integer,intent(in)           :: nS
  double precision,intent(in)  :: ENuc
  complex*16,intent(in)        :: ERHF
  complex*16,intent(in)        :: eHF(nBas)
  complex*16,intent(in)        :: ERI(nBas,nBas,nBas,nBas)
  complex*16,intent(in)        :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                      :: ia
  integer                      :: ispin
  logical                      :: dRPA
  double precision             :: lambda
  complex*16,allocatable       :: Aph(:,:)
  complex*16,allocatable       :: Bph(:,:)
  complex*16,allocatable       :: Om(:)
  complex*16,allocatable       :: XpY(:,:)
  complex*16,allocatable       :: XmY(:,:)

  complex*16                   :: EcRPA(nspin)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted ph-RPA Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA = .true.
  EcRPA(:) = 0d0
  lambda = 1d0

! Memory allocation

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS),Aph(nS,nS),Bph(nS,nS))

! Singlet manifold

  if(singlet) then 

    ispin = 1

                 call complex_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,Aph)
    if(.not.TDA) call complex_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

    call complex_phRLR(TDA,nS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPA@RHF','singlet',nS,Om)
    !call complex_phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  end if

! Triplet manifold 

  if(triplet) then 

    ispin = 2

                 call complex_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,Aph)
    if(.not.TDA) call complex_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

    call complex_phRLR(TDA,nS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPA@RHF','triplet',nS,Om)
    !call complex_phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if
  
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPA correlation energy (spin-conserved) = ',real(EcRPA(1)),'+ i*',aimag(EcRPA(1)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPA correlation energy (spin-flip)      = ',real(EcRPA(2)),'+ i*',aimag(EcRPA(2)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPA correlation energy                  = ',real(sum(EcRPA)),'+ i*',aimag(sum(EcRPA)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPA total energy                        = ',real(ENuc + ERHF + sum(EcRPA)),'+ i*'&
                                                                                            ,aimag(ENuc + ERHF + sum(EcRPA)),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  deallocate(Om,XpY,XmY,Aph,Bph)

end subroutine 
