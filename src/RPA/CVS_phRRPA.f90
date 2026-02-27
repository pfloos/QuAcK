subroutine CVS_phRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,nCVS,FC,occupations,ENuc,ERHF,ERI,dipole_int,eHF)

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
  integer,intent(in)           :: nCVS
  integer,intent(in)           :: FC
  integer,intent(in)           :: occupations(nO)
  double precision,intent(in)  :: ENuc
  double precision,intent(in)  :: ERHF
  double precision,intent(in)  :: eHF(nBas)
  double precision,intent(in)  :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)  :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                      :: ia,i
  integer                      :: ispin
  integer                      :: nSCVS, nFC
  logical                      :: dRPA,found
  double precision             :: t1, t2
  double precision             :: lambda
  double precision,allocatable :: Aph(:,:)
  double precision,allocatable :: Bph(:,:)
  double precision,allocatable :: Om(:)
  double precision,allocatable :: XpY(:,:)
  double precision,allocatable :: XmY(:,:)
  integer,allocatable          :: virtuals(:)
  integer,allocatable          :: occupations_fc(:)

  double precision             :: EcRPA(nspin)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted ph-RPA Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! CVS

  print *, "No exications to the first", nCVS, "orbital(s) are considered."
  if(nC/=0) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core

  nFC = MERGE(1,0,FC/=0) 
  allocate(occupations_fc(nO-nFC))
  ! remove FC from occupations
  do ispin=1,nspin
    occupations_fc(1:nO-nFC) = occupations(1:nO - nFC) 
    found = .false.
    do i=1,nO-1
      if(.not. found) then
        if(occupations(i)==FC) then
          found = .true.
          occupations_fc(i) = occupations(i+1) 
        else
          occupations_fc(i) = occupations(i)
        endif
      else
        occupations_fc(i) = occupations(i+1) 
      endif 
    enddo
  enddo
  print *, "Not Frozen orbitals:"
  print *,occupations_fc(1:nO-nFC)

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
  if(nC/=0) then
    print *, "No frozen core with CVS available !"
    stop
  end if
  nSCVS = (nV - nCVS)*(nO - nFC)
  allocate(Om(nSCVS),XpY(nSCVS,nSCVS),XmY(nSCVS,nSCVS),Aph(nSCVS,nSCVS),Bph(nSCVS,nSCVS),virtuals(nV))
  call non_occupied(nO,nBas,occupations,virtuals)

! Singlet manifold

  if(singlet) then 

    ispin = 1

                 call CVS_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,eHF,ERI,Aph)
    if(.not.TDA) call CVS_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,ERI,Bph)

    call CVS_phRLR(TDA,nSCVS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('phRPA@RHF','singlet',nSCVS,Om)
    call CVS_phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,dipole_int,Om,XpY,XmY)

  end if

! Triplet manifold 

  if(triplet) then 

    ispin = 2

                 call CVS_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,eHF,ERI,Aph)
    if(.not.TDA) call CVS_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,ERI,Bph)

    call CVS_phRLR(TDA,nSCVS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('phRPA@RHF','triplet',nSCVS,Om)
    call CVS_phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,dipole_int,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy (singlet) = ',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy (triplet) = ',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy           = ',sum(EcRPA),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF total energy                 = ',ENuc + ERHF + sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  deallocate(Om,XpY,XmY,Aph,Bph)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then
    print *,"No ADC implemented with CVS."
  end if

  if(dotest) then

    call dump_test_value('R','phRPA correlation energy',sum(EcRPA))

  end if

end subroutine 
