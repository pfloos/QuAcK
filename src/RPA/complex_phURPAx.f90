subroutine complex_phURPAx(dotest,TDA,doACFDT,exchange_kernel,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF, & 
                  ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,CAP_MO,eHF,c,S,occupations)

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
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  complex*16,intent(in)         :: EUHF
  complex*16,intent(in)         :: eHF(nBas,nspin)
  complex*16,intent(in)         :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  complex*16,intent(in)         :: ERI_aaaa(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: ERI_aabb(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: ERI_bbbb(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: dipole_int_aa(nBas,nBas,ncart)
  complex*16,intent(in)         :: dipole_int_bb(nBas,nBas,ncart)
  complex*16,intent(in)         :: CAP_MO(nBas,nBas,nspin)

! Local variables

  logical                       :: dRPA,found
  integer                       :: ispin,i
  double precision              :: lambda

  integer                       :: nSa,nSb,nSt,nFC(nspin)
  complex*16,allocatable        :: Aph(:,:)
  complex*16,allocatable        :: Bph(:,:)
  complex*16,allocatable        :: Om(:)
  complex*16,allocatable        :: XpY(:,:)
  complex*16,allocatable        :: XmY(:,:)
  integer,allocatable           :: virtuals(:,:)
  integer,allocatable           :: occupations_fc(:,:)


  complex*16                    :: EcRPA(nspin)

! Hello world
 
  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Unrestricted ph-RPA Calculation *'
  write(*,*)'***********************************'
  write(*,*)

! CVS
  if(any(nCVS/=0)) then
  print *, "CVS is applied."
    print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
    print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
    if(any(nC/=0)) then
      print *, "Do not use PyDuck frozen core with CVS !"
      stop
    end if
  end if

! Frozen Core

  nFC(1) = MERGE(1,0,FC(1)/=0) 
  nFC(2) = MERGE(1,0,FC(2)/=0)
  
  allocate(occupations_fc(maxval(nO-nFC),nspin))
  if(any(nFC /= 0)) then 
    ! remove FC from occupations
    do ispin=1,nspin
      occupations_fc(1:nO(ispin)-nFC(ispin),ispin) = occupations(1:nO(ispin) - nFC(ispin), ispin) 
      found = .false.
      do i=1,nO(ispin)-1
        if(.not. found) then
          if(occupations(i,ispin)==FC(ispin)) then
            found = .true.
            occupations_fc(i,ispin) = occupations(i+1,ispin) 
          else
            occupations_fc(i,ispin) = occupations(i,ispin)
          endif
        else
          occupations_fc(i,ispin) = occupations(i+1,ispin) 
        endif 
      enddo
      print *, "Frozen Core is applied." 
      print *, "Not Frozen orbitals:"
      print *,occupations_fc(1:nO(ispin)-nFC(ispin),ispin)
    enddo
  else
    do ispin=1,nspin
     occupations_fc(1:nO(ispin),ispin) = occupations(1:nO(ispin),ispin) 
    end do
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
  allocate(virtuals(nBas - minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do

! Spin-conserved transitions

  if(spin_conserved) then 

    ispin = 1

    ! Memory allocation
    nFC(1) = MERGE(1,0,FC(1)/=0) 
    nFC(2) = MERGE(1,0,FC(2)/=0) 
    nSa = (nBas - nO(1) - nCVS(1))*(nO(1) - nFC(1))
    nSb = (nBas - nO(2) - nCVS(2))*(nO(2) - nFC(2))
    nSt = nSa + nSb

    allocate(Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt))

    call complex_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc, virtuals,&
                            lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call complex_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc, virtuals,&
                             lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

    call complex_phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPAx@UHF','spin-conserved',nSt,Om)
    call complex_phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,&
                nCVS,nFC,occupations_fc,virtuals,dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

    deallocate(Aph,Bph,Om,XpY,XmY)

  end if

! Spin-flip transitions

  if(spin_flip) then

    ispin = 2
    
    ! Memory allocation
    
    nSa = (nBas - nO(2) - nCVS(2))*(nO(1) - nFC(1))
    nSb = (nBas - nO(1) - nCVS(1))*(nO(2) - nFC(2))
    nSt = nSa + nSb

    allocate(Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt))

    call complex_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,&
                     lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call complex_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,&
                                  lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

    call complex_phULR(TDA,nSa,nSa,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call complex_print_excitation_energies('phRPAx@UHF','spin-flip',nSt,Om)
    call complex_phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,&
                nCVS,nFC,occupations_fc,virtuals,dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

    deallocate(Aph,Bph,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 0.5d0*EcRPA(2)

  else

    EcRPA(2) = 0d0

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phURPAx correlation energy (spin-conserved) = ',real(EcRPA(1)),'+ i*',aimag(EcRPA(1)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phURPAx correlation energy (spin-flip)      = ',real(EcRPA(2)),'+ i*',aimag(EcRPA(2)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phURPAx correlation energy                  = ',real(sum(EcRPA)),'+ i*',aimag(sum(EcRPA)),' au'
  write(*,'(2X,A50,F20.10,A3,F20.10,A3)') 'Tr@phURPAx total energy                        = ',real(ENuc + EUHF + sum(EcRPA)),'+ i*',aimag(ENuc + EUHF + sum(EcRPA)),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then
    print *, "No adiabatic connections yet with CVS !"
  end if

end subroutine 
