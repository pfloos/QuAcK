subroutine complex_phRRPAx(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,nCVS,FC,occupations,ENuc,ERHF,ERI,dipole_int,CAP_MO,eHF)

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
  complex*16,intent(in)        :: ERHF
  complex*16,intent(in)        :: eHF(nBas)
  complex*16,intent(in)        :: ERI(nBas,nBas,nBas,nBas)
  complex*16,intent(in)        :: dipole_int(nBas,nBas,ncart)
  complex*16,intent(in)        :: CAP_MO(nBas,nBas)

! Local variables

  integer                      :: ia,nSCVS,nFC,i
  integer                      :: ispin
  integer,allocatable          :: virtuals(:)
  integer,allocatable          :: occupations_fc(:)
  logical                      :: dRPA, found
  double precision             :: lambda
  complex*16,allocatable       :: Aph(:,:)
  complex*16,allocatable       :: Bph(:,:)
  complex*16,allocatable       :: Om(:)
  complex*16,allocatable       :: XpY(:,:)
  complex*16,allocatable       :: XmY(:,:)

  complex*16                   :: EcRPA(nspin)
  complex*16                   :: EcW(nspin)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted ph-RPA Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! CVS
  if(nCVS/=0) then
  print *, "CVS is applied."
    print *, "No exications to the first", nCVS, " spatial orbital(s) are considered."
    if(nC/=0) then
      print *, "Do not use PyDuck frozen core with CVS !"
      stop
    end if
  end if

! Frozen Core

  nFC = MERGE(1,0,FC/=0) 
  allocate(occupations_fc(nO-nFC))
  if(nFC/=0) then
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
  else
    occupations_fc(1:nO) = occupations(1:nO)
  end if

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA = .false.
  EcRPA(:) = 0d0
  lambda = 1d0
  nSCVS = (nV - nCVS)*(nO-nFC)

! Memory allocation

  allocate(Om(nSCVS),XpY(nSCVS,nSCVS),XmY(nSCVS,nSCVS),Aph(nSCVS,nSCVS),Bph(nSCVS,nSCVS),virtuals(nV))
  call non_occupied(nO,nBas,occupations(1:nO),virtuals)

! Singlet manifold

  if(singlet) then 

    ispin = 1

                 call complex_CVS_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,eHF,ERI,Aph)
    if(.not.TDA) call complex_CVS_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,ERI,Bph)

    call complex_phRLR(TDA,nSCVS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPAx@RHF','singlet',nS,Om)
    call complex_phLR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,&
                nCVS,nFC,occupations_fc,virtuals,Om,XpY,XmY)

  end if

! Triplet manifold 

  if(triplet) then 

    ispin = 2

                 call complex_CVS_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,eHF,ERI,Aph)
    if(.not.TDA) call complex_CVS_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,lambda,ERI,Bph)

    call complex_phRLR(TDA,nSCVS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPAx@RHF','triplet',nSCVS,Om)
    call complex_phLR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,&
                nCVS,nFC,occupations_fc,virtuals,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if
  
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPAx correlation energy (spin-conserved) = ',real(EcRPA(1)),'+ i*',aimag(EcRPA(1)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPAx correlation energy (spin-flip)      = ',real(EcRPA(2)),'+ i*',aimag(EcRPA(2)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPAx correlation energy                  = ',real(sum(EcRPA)),'+ i*',aimag(sum(EcRPA)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phRRPAx total energy                        = ',real(ENuc + ERHF + sum(EcRPA)),'+ i*'&
                                                                                            ,aimag(ENuc + ERHF + sum(EcRPA)),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  deallocate(Om,XpY,XmY,Aph,Bph,virtuals)

end subroutine 
