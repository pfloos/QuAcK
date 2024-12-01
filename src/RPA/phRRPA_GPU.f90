#ifdef USE_GPU

subroutine phRRPA_GPU(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

  use cu_quack_module


  implicit none
  include 'parameters.h'
  include 'quadrature.h'


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
  double precision,intent(in)  :: ERHF
  double precision,intent(in)  :: eHF(nBas)
  double precision,intent(in)  :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)  :: dipole_int(nBas,nBas,ncart)


  integer                      :: i, a, ia
  integer                      :: ispin
  logical                      :: dRPA
  double precision             :: t1, t2
  integer,         allocatable :: iorder(:)
  double precision,allocatable :: Om(:)
  double precision,allocatable :: XpY(:,:)
  double precision,allocatable :: XmY(:,:)

  double precision             :: EcRPA(nspin)


  write(*,*)
  write(*,*)'******************************************'
  write(*,*)'* Restricted ph-RPA Calculation (on GPU) *'
  write(*,*)'******************************************'
  write(*,*)

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA = .true.
  EcRPA(:) = 0d0


  allocate(Om(nS), XpY(nS,nS), XmY(nS,nS))

  if(singlet) then 

    if(TDA) then

      !print*, 'start diag on GPU:'
      !call wall_time(t1)
      call ph_drpa_tda_sing(nO, nBas, nS, eHF(1), ERI(1,1,1,1), Om(1), XpY(1,1))
      !call wall_time(t2)
      !print*, 'diag time on GPU (sec):', t2 - t1
      XmY(:,:) = XpY(:,:)

    else

      !call wall_time(t1)
      call ph_drpa_sing(nO, nBas, nS, eHF(1), ERI(1,1,1,1), Om(1), XpY(1,1), XmY(1,1))
      call wall_time(t2)
      print *, "wall time for dRPA on GPU (sec) = ", t2 - t1
      !do ia = 1, nS
      !  write(111, *) Om(ia)
      !enddo
      !stop

    endif

    ! TODO
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = transpose(XmY(:,:))

    call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)
  endif

  if(triplet) then 

    XpY(:,:) = 0.d0
    allocate(iorder(nS))
    ia = 0
    do i = nC+1, nO
      do a = nO+1, nBas-nR
        ia = ia + 1
        iorder(ia) = ia
        Om(ia) = eHF(a) - eHF(i)
        XpY(ia,ia) = 1.d0
      enddo
    enddo

    call quick_sort(Om(1), iorder(1), nS)
    deallocate(iorder)

    XmY(:,:) = XpY(:,:)

    call print_excitation_energies('phRPA@RHF','triplet',nS,Om)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)
  endif

  deallocate(Om, XpY, XmY)


  ! TODO
  ! init EcRPA
  if(exchange_kernel) then
    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)
  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy (singlet) = ',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy (triplet) = ',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF correlation energy           = ',sum(EcRPA),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPA@RHF total energy                 = ',ENuc + ERHF + sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    ! TODO
    call phACFDT(exchange_kernel,dRPA,TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,eHF,EcRPA)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'AC@phRPA@RHF correlation energy (singlet) = ',EcRPA(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@phRPA@RHF correlation energy (triplet) = ',EcRPA(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@phRPA@RHF correlation energy           = ',sum(EcRPA),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@phRPA@RHF total energy                 = ',ENuc + ERHF + sum(EcRPA),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  endif

  if(dotest) then
    call dump_test_value('R','phRPA correlation energy (on GPU)',sum(EcRPA))
  endif

end subroutine 

#else

subroutine phRRPA_GPU(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

  logical,intent(in)            :: dotest
  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  print*, "compile with USE_GPU FLAG!"
  stop
end

#endif
