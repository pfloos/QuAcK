subroutine ppURPA(TDA,doACFDT,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,e)

! Perform unrestricted pp-RPA calculations

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,allocatable  :: Omega1sc(:),Omega1sf(:)
  double precision,allocatable  :: X1sc(:,:),X1sf(:,:)
  double precision,allocatable  :: Y1sc(:,:),Y1sf(:,:)
  double precision,allocatable  :: Omega2sc(:),Omega2sf(:)
  double precision,allocatable  :: X2sc(:,:),X2sf(:,:)
  double precision,allocatable  :: Y2sc(:,:),Y2sf(:,:)

  double precision              :: Ec_ppURPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'|  particle-particle URPA calculation  |'
  write(*,*)'****************************************'
  write(*,*)

! Initialization

  Ec_ppURPA(:) = 0d0
  EcAC(:)   = 0d0

! Spin-conserved manifold

  if(spin_conserved) then 

    ispin = 1
    iblock = 1

!Spin-conserved quantities

    nPab = nV(1)*nV(2)
    nHab = nO(1)*nO(2)

    nP_sc = nPab
    nH_sc = nHab

! Memory allocation

    allocate(Omega1sc(nP_sc),X1sc(nP_sc,nP_sc),Y1sc(nH_sc,nP_sc), &
             Omega2sc(nH_sc),X2sc(nP_sc,nH_sc),Y2sc(nH_sc,nH_sc))

    call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                         nP_sc,nHaa,nHab,nHbb,nH_sc,1d0,e,ERI_aaaa, &
                                         ERI_aabb,ERI_bbbb,Omega1sc,X1sc,Y1sc, &
                                         Omega2sc,X2sc,Y2sc,Ec_ppURPA(ispin))
    
    call print_excitation('pp-RPA (N+2)',5,nP_sc,Omega1sc)
    call print_excitation('pp-RPA (N-2)',5,nH_sc,Omega2sc)

  endif

! Spin-flip manifold 

  if(spin_flip) then 

    ispin = 2
    iblock = 2 

!Spin-flip quantities

    nPaa = nV(1)*(nV(1)-1)/2
    nPbb = nV(2)*(nV(2)-1)/2

    nP_sf  = nPaa 

    nHaa = nO(1)*(nO(1)-1)/2
    nHbb = nO(2)*(nO(2)-1)/2

    nH_sf  = nHaa 

    allocate(Omega1sf(nP_sf),X1sf(nP_sf,nP_sf),Y1sf(nH_sf,nP_sf), &
             Omega2sf(nH_sf),X2sf(nP_sf,nH_sf),Y2sf(nH_sf,nH_sf))

    call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                         nP_sf,nHaa,nHab,nHbb,nH_sf,1d0,e,ERI_aaaa, &
                                         ERI_aabb,ERI_bbbb,Omega1sf,X1sf,Y1sf, &
                                         Omega2sf,X2sf,Y2sf,Ec_ppURPA(ispin))

    deallocate(Omega1sf,X1sf,Y1sf,Omega2sf,X2sf,Y2sf)

    iblock = 3

    nP_sf  = nPbb
    nH_sf  = nHbb

    allocate(Omega1sf(nP_sf),X1sf(nP_sf,nP_sf),Y1sf(nH_sf,nP_sf), &
             Omega2sf(nH_sf),X2sf(nP_sf,nH_sf),Y2sf(nH_sf,nH_sf))

    call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,&
                                         nP_sf,nHaa,nHab,nHbb,nH_sf,1d0,e,ERI_aaaa,&
                                         ERI_aabb,ERI_bbbb,Omega1sf,X1sf,Y1sf,&
                                         Omega2sf,X2sf,Y2sf,Ec_ppURPA(ispin))

    call print_excitation('pp-RPA (N+2)',6,nP_sf,Omega1sf)
    call print_excitation('pp-RPA (N-2)',6,nH_sf,Omega2sf)

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (spin-conserved) =',Ec_ppURPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (spin-flip)      =',3d0*Ec_ppURPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy                  =',Ec_ppURPA(1) + 3d0*Ec_ppURPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA total energy                        =',ENuc + EUHF + Ec_ppURPA(1) + 3d0*Ec_ppURPA(2)
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


end subroutine ppURPA
