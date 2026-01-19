subroutine CVS_phURPA(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,nCVS,ENuc,EUHF, & 
                  ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF,c,S,occupations)

! Perform random phase approximation calculation with exchange (aka TDHF) in the unrestricted formalism

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: lambda

  integer                       :: nSa,nSb,nSt
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  integer,allocatable           :: virtuals(:,:)


  double precision              :: rho
  double precision              :: EcRPA(nspin)

! Hello world
 
  write(*,*)
  write(*,*)'***************************************'
  write(*,*)'* Unrestricted CVS ph-RPA Calculation *'
  write(*,*)'***************************************'
  write(*,*)

! CVS

  print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
  print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
  if(any(nC/=0)) then
    print *, "Do not use frozen core with CVS !"
    stop
  end if

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA = .true.
  EcRPA(:) = 0d0
  lambda = 1d0
  allocate(virtuals(nBas - minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do

! Spin-conserved transitions

  if(spin_conserved) then 

    ispin = 1

    ! Memory allocation

    nSa = (nBas - nO(1) - nCVS(1))*nO(1)
    nSb = (nBas - nO(2) - nCVS(2))*nO(2) 
    nSt = nSa + nSb

    allocate(Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt))

    call CVS_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS, occupations, virtuals,&
                            lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call CVS_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS, occupations, virtuals,&
                             lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

    call phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('phRPA@UHF','spin-conserved',nSt,Om)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt, &
                               dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

    deallocate(Aph,Bph,Om,XpY,XmY)

  end if

! Spin-flip transitions

  if(spin_flip) then

    ispin = 2
    print *, "Spin flip transitions not yet implemented for CVS !"
    ! Memory allocation
    
    nSa = (nBas - nO(2) - nCVS(2))*nO(1)
    nSb = (nBas - nO(1) - nCVS(1))*nO(2) 
    nSt = nSa + nSb

    allocate(Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt))

    call phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

    call phULR(TDA,nSa,nSa,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('phRPA@UHF','spin-flip',nSt,Om)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

    deallocate(Aph,Bph,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phURPA correlation energy (spin-conserved) = ',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phURPA correlation energy (spin-flip)      = ',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phURPA correlation energy                  = ',sum(EcRPA),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phURPA total energy                        = ',ENuc + EUHF + sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then
    print *, "No adiabatic connections yet with CVS !"
  end if

  if(dotest) then

    call dump_test_value('U','phRPA correlation energy',sum(EcRPA))

  end if

end subroutine 
