subroutine ppURPA(dotest,TDA,doACFDT,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,e)

! Perform unrestricted pp-RPA calculations

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin,iblock
  integer                       :: nPaa,nPbb,nPab,nP_sc,nP_sf
  integer                       :: nHaa,nHbb,nHab,nH_sc,nH_sf
  double precision,allocatable  :: Om1sc(:),Om1sf(:)
  double precision,allocatable  :: X1sc(:,:),X1sf(:,:)
  double precision,allocatable  :: Y1sc(:,:),Y1sf(:,:)
  double precision,allocatable  :: Om2sc(:),Om2sf(:)
  double precision,allocatable  :: X2sc(:,:),X2sf(:,:)
  double precision,allocatable  :: Y2sc(:,:),Y2sf(:,:)

  double precision              :: EcRPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Unrestricted pp-RPA Calculation *'
  write(*,*)'***********************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0
  EcAC(:)   = 0d0

!alpha-beta block  

  ispin = 1
  iblock = 3

  nPab = nV(1)*nV(2)
  nHab = nO(1)*nO(2)

  nP_sc = nPab
  nH_sc = nHab

! Memory allocation

  allocate(Om1sc(nP_sc),X1sc(nP_sc,nP_sc),Y1sc(nH_sc,nP_sc), &
           Om2sc(nH_sc),X2sc(nP_sc,nH_sc),Y2sc(nH_sc,nH_sc))

  call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
             nP_sc,nHaa,nHab,nHbb,nH_sc,1d0,e,ERI_aaaa, &
             ERI_aabb,ERI_bbbb,Om1sc,X1sc,Y1sc, &
             Om2sc,X2sc,Y2sc,EcRPA(ispin))
    
  call print_excitation_energies('ppRPA@UHF','2p (alpha-beta)',nP_sc,Om1sc)
  call print_excitation_energies('ppRPA@UHF','2h (alpha-beta)',nH_sc,Om2sc)

!alpha-alpha block

  ispin = 2
  iblock = 4

  nPaa = nV(1)*(nV(1)-1)/2
  nPbb = nV(2)*(nV(2)-1)/2

  nP_sf  = nPaa 

  nHaa = nO(1)*(nO(1)-1)/2
  nHbb = nO(2)*(nO(2)-1)/2

  nH_sf  = nHaa 

  allocate(Om1sf(nP_sf),X1sf(nP_sf,nP_sf),Y1sf(nH_sf,nP_sf), &
           Om2sf(nH_sf),X2sf(nP_sf,nH_sf),Y2sf(nH_sf,nH_sf))

  call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
             nP_sf,nHaa,nHab,nHbb,nH_sf,1d0,e,ERI_aaaa, &
             ERI_aabb,ERI_bbbb,Om1sf,X1sf,Y1sf, &
             Om2sf,X2sf,Y2sf,EcRPA(ispin))

  call print_excitation_energies('ppRPA@UHF','2h (alpha-alpha)',nP_sf,Om1sf)
  call print_excitation_energies('ppRPA@UHF','2p (alpha-alpha)',nH_sf,Om2sf)

  deallocate(Om1sf,X1sf,Y1sf,Om2sf,X2sf,Y2sf)

!beta-beta block

  iblock = 7 

  nP_sf  = nPbb
  nH_sf  = nHbb

  allocate(Om1sf(nP_sf),X1sf(nP_sf,nP_sf),Y1sf(nH_sf,nP_sf), &
           Om2sf(nH_sf),X2sf(nP_sf,nH_sf),Y2sf(nH_sf,nH_sf))

  call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,&
             nP_sf,nHaa,nHab,nHbb,nH_sf,1d0,e,ERI_aaaa,&
             ERI_aabb,ERI_bbbb,Om1sf,X1sf,Y1sf,&
             Om2sf,X2sf,Y2sf,EcRPA(ispin))

  call print_excitation_energies('ppRPA@UHF','2p (beta-beta)',nP_sf,Om1sf)
  call print_excitation_energies('ppRPA@UHF','2h (beta-beta)',nH_sf,Om2sf)

  EcRPA(2) = 3d0*EcRPA(2)

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppURPA correlation energy (spin-conserved) = ',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppURPA correlation energy (spin-flip)      = ',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppURPA correlation energy                  = ',sum(EcRPA),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppURPA total energy                        = ',ENuc + EUHF + sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

! if(doACFDT) then

!   write(*,*) '---------------------------------------------------------'
!   write(*,*) 'Adiabatic connection version of pp-RPA correlation energy'
!   write(*,*) '---------------------------------------------------------'
!   write(*,*)

!   call ACFDT_pp(TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,e,EcAC)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (singlet) =',EcAC(1),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (triplet) =',EcAC(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy           =',EcAC(1) + EcAC(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2),' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('U','ppRPA correlation energy',sum(EcRPA))

  end if

end subroutine 
