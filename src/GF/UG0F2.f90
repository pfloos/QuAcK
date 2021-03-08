subroutine UG0F2(BSE,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,linearize,eta,nBas,nC,nO,nV,nR,nS,ENuc,EUHF, &
                 S,ERI,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF)

! Perform unrestricted G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas,nspin)

! Local variables

  integer                       :: is
  integer                       :: ispin
  double precision              :: Ec(nsp)
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc

  double precision,allocatable  :: eGF2lin(:,:)
  double precision,allocatable  :: eGF2(:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0F2 calculation           |'
  write(*,*)'|         *** Unrestricted version ***         |'
  write(*,*)'************************************************'
  write(*,*)

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),eGF2(nBas,nspin),eGF2lin(nBas,nspin))

!---------------------!
! Compute self-energy !
!---------------------!

  call unrestricted_self_energy_GF2_diag(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGF2,SigC,Z)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGF2lin(:,:) = eHF(:,:) + Z(:,:)*SigC(:,:)

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGF2(:,:) = eGF2lin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

    print*,'!!! Graphical solution NYI for UG0F2 !!!'
    stop
 
  end if

! Compute MP2 correlation energy

  call UMP2(nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EUHF,eHF,Ec)

! Dump results

  call print_UG0F2(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGF2,Ec)

! Perform BSE calculation

  if(BSE) then

    print*,'!!! BSE2 NYI for UG0F2 !!!'

  end if

end subroutine UG0F2
