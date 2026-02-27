subroutine CVS_UG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
                 linearize,eta,doSRG,nBas,nC,nO,nV,nR,nS,nCVS,FC,ENuc,EUHF,S,ERI_aaaa,ERI_aabb,ERI_bbbb,       & 
                 dipole_int_aa,dipole_int_bb,cHF,eHF,occupations)

! Perform unrestricted G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)

! Local variables

  logical                       :: print_W = .true.
  logical                       :: dRPA_W 
  logical                       :: found
  integer                       :: is,i,ispin
  integer                       :: isp_W
  double precision              :: flow
  double precision              :: EcRPA
  double precision              :: EcGM(nspin)
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nSa,nSb,nSt,nFC(nspin)
  integer,allocatable           :: virtuals(:,:)
  integer,allocatable           :: occupations_fc(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:,:)

  double precision,allocatable  :: eGWlin(:,:)
  double precision,allocatable  :: eGW(:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Unrestricted G0W0 Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! CVS

  print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
  print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
  if(any(nC/=0)) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core

  nFC(1) = MERGE(1,0,FC(1)/=0) 
  nFC(2) = MERGE(1,0,FC(2)/=0)
  allocate(occupations_fc(maxval(nO-nFC),nspin))
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
  enddo
  do ispin=1,nspin
    print *, "Not Frozen orbitals:"
    print *,occupations_fc(1:nO(ispin)-nFC(ispin),ispin)
  end do

! Initialization

  EcRPA = 0d0
  dRPA_W = .true.
  
  allocate(virtuals(nBas - minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G0W0 scheme ***'
    write(*,*)

  end if


! Memory allocation

  nSa = (nBas - nO(1) - nCVS(1))*(nO(1) - nFC(1))
  nSb = (nBas - nO(2) - nCVS(2))*(nO(2) - nFC(2))
  nSt = nSa + nSb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),eGWlin(nBas,nspin),eGW(nBas,nspin), &
           Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt),rho(nBas,nBas,nSt,nspin))

!-------------------!
! Compute screening !
!-------------------!

! Spin-conserving transitions

  isp_W = 1

                 call CVS_phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
  if(.not.TDA_W) call CVS_phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

  call CVS_phULR(TDA_W,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)
  
  if(print_W) call print_excitation_energies('phRPA@UHF','spin-conserved',nSt,Om)

!----------------------!
! Excitation densities !
!----------------------!
 
  call CVS_UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

!------------------------------------------------!
! Compute self-energy and renormalization factor !
!------------------------------------------------!

  if(doSRG) then
    call UGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nSt,eHF,Om,rho,EcGM,SigC,Z)
  else
    call UGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nSt,eHF,Om,rho,EcGM,SigC,Z)
  end if

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGWlin(:,:) = eHF(:,:) + Z(:,:)*SigC(:,:) 

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:,:) = eGWlin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)

    do is=1,nspin

      write(*,*)'-----------------------------------------------------'
      if(is==1) write(*,*)'    Spin-up   orbitals    '
      if(is==2) write(*,*)'    Spin-down orbitals    '

      call CVS_UGW_QP_graph(doSRG,eta,flow,nBas,nC(is),nO(is),nV(is),nR(is),nSt,nCVS(is),nFC(is),occupations_fc(:,is),virtuals(:,is),eHF(:,is), & 
                        Om,rho(:,:,:,is),eGWlin(:,is),eHF(:,is),eGW(:,is),Z(:,is))
    end do
 
  end if

! Compute RPA correlation energy

                 call CVS_phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
  if(.not.TDA_W) call CVS_phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
    
  call CVS_phULR(TDA_W,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

! Dump results

  call print_UG0W0(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA,EcGM)

! Free memory

  deallocate(Om,XpY,XmY,rho)

! Perform BSE calculation
!
!  if(dophBSE) then
!
!    call UGW_phBSE(exchange_kernel,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS,S, &
!                   ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eHF,eGW,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A60,F20.10,A3)') 'Tr@BSE@G0W0@UHF correlation energy (spin-conserved) = ',EcBSE(1),' au'
!    write(*,'(2X,A60,F20.10,A3)') 'Tr@BSE@G0W0@UHF correlation energy (spin-flip)      = ',EcBSE(2),' au'
!    write(*,'(2X,A60,F20.10,A3)') 'Tr@BSE@G0W0@UHF correlation energy                  = ',sum(EcBSE),' au'
!    write(*,'(2X,A60,F20.10,A3)') 'Tr@BSE@G0W0@UHF total       energy                  = ',ENuc + EUHF + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!!   Compute the BSE correlation energy via the adiabatic connection 
!
!    if(doACFDT) then
!
!      call UGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,spin_conserved,spin_flip,eta, & 
!                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGW,EcBSE)
!
!      write(*,*)
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,'(2X,A60,F20.10,A3)') 'AC@BSE@G0W0@UHF correlation energy (spin-conserved) = ',EcBSE(1),' au'
!      write(*,'(2X,A60,F20.10,A3)') 'AC@BSE@G0W0@UHF correlation energy (spin-flip)      = ',EcBSE(2),' au'
!      write(*,'(2X,A60,F20.10,A3)') 'AC@BSE@G0W0@UHF correlation energy                  = ',sum(EcBSE),' au'
!      write(*,'(2X,A60,F20.10,A3)') 'AC@BSE@G0W0@UHF total       energy                  = ',ENuc + EUHF + sum(EcBSE),' au'
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,*)
!
!    end if
!
!  end if
!
!! Testing zone
!
!  if(dotest) then
!
!    call dump_test_value('U','G0W0 correlation energy',EcRPA)
!    call dump_test_value('U','G0W0 HOMOa energy',eGW(nO(1),1))
!    call dump_test_value('U','G0W0 LUMOa energy',eGW(nO(1)+1,1))
!    call dump_test_value('U','G0W0 HOMOa energy',eGW(nO(2),2))
!    call dump_test_value('U','G0W0 LUMOa energy',eGW(nO(2)+1,2))
!
!  end if

end subroutine 
