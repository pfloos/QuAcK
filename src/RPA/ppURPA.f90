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

  integer                       :: ispin
  integer                       :: nS
  integer                       :: nPaa,nPbb,nPab,nPt
  integer                       :: nHaa,nHbb,nHab,nHt
  double precision,allocatable  :: Omega1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: Omega2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

  double precision              :: Ec_ppRPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'|  particle-particle URPA calculation  |'
  write(*,*)'****************************************'
  write(*,*)

! Initialization

  Ec_ppRPA(:) = 0d0
  EcAC(:)   = 0d0

! Useful quantities

  nPaa = nV(1)*(nV(1)-1)/2
  nPab = nV(1)*nV(2)
  nPbb = nV(2)*nV(2)
  nPt  = nPaa + nPab + nPbb
 
  nHaa = nO(1)*(nO(1)-1)/2
  nHab = nO(1)*nO(2)
  nHbb = nO(2)*nO(2)
  nHt  = nHaa + nHab + nHbb

 ! Memory allocation

 allocate(Omega1(nPt),X1(nPt,nPt),Y1(nHt,nPt), & 
          Omega2(nHt),X2(nPt,nHt),Y2(nHt,nHt))

! allocate(Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
!          Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt))

! Spin-conserved manifold

  if(spin_conserved) then 

    ispin = 1

!   call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,e,ERI, & 
!                           Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,Ec_ppRPA(ispin))

!   call print_excitation('pp-RPA (N+2)',5,nVVs,Omega1s)
!   call print_excitation('pp-RPA (N-2)',5,nOOs,Omega2s)

  endif

! Spin-flip manifold 

  if(spin_flip) then 

    ispin = 2

!   call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,e,ERI, &
!                           Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,Ec_ppRPA(ispin))

!   call print_excitation('pp-RPA (N+2)',6,nVVt,Omega1t)
!   call print_excitation('pp-RPA (N-2)',6,nOOt,Omega2t)

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (spin-conserved) =',Ec_ppRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (spin-flip)      =',3d0*Ec_ppRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy                  =',Ec_ppRPA(1) + 3d0*Ec_ppRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA total energy                        =',ENuc + EUHF + Ec_ppRPA(1) + 3d0*Ec_ppRPA(2)
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
