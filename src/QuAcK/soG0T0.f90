subroutine soG0T0(eta,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,eHF)

! Perform G0W0 calculation with a T-matrix self-energy (G0T0) in the spinorbital basis

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: nOO
  integer                       :: nVV
  double precision              :: EcRPA
  integer                       :: nBas2,nC2,nO2,nV2,nR2
  double precision,allocatable  :: Omega1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: rho1(:,:,:)
  double precision,allocatable  :: Omega2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)
  double precision,allocatable  :: rho2(:,:,:)
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: eG0T0(:)
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)

! Output variables


! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot soG0T0 calculation         |'
  write(*,*)'************************************************'
  write(*,*)

! Spatial to spin orbitals

  nBas2 = 2*nBas

  allocate(seHF(nBas2),sERI(nBas2,nBas2,nBas2,nBas2))

  call spatial_to_spin_MO_energy(nBas,eHF,nBas2,seHF)
  call spatial_to_spin_ERI(nBas,ERI,nBas2,sERI)

! Define occupied and virtual spaces

  nO2 = 2*nO
  nV2 = 2*nV
  nC2 = 2*nC 
  nR2 = 2*nR

! Dimensions of the rr-RPA linear reponse matrices

  nOO = nO2*(nO2-1)/2
  nVV = nV2*(nV2-1)/2

! Memory allocation

  allocate(Omega1(nVV),X1(nVV,nVV),Y1(nOO,nVV), & 
           Omega2(nOO),X2(nVV,nOO),Y2(nOO,nOO), & 
           rho1(nBas2,nBas2,nVV),rho2(nBas2,nBas2,nOO), & 
           eG0T0(nBas2),SigT(nBas2),Z(nBas2))

!----------------------------------------------
! Spinorbital basis
!----------------------------------------------

 ispin = 2

! Compute linear response

  call linear_response_pp(ispin,.false.,nBas2,nC2,nO2,nV2,nR2, & 
                          nOO,nVV,seHF(:),sERI(:,:,:,:),  & 
                          Omega1(:),X1(:,:),Y1(:,:),   & 
                          Omega2(:),X2(:,:),Y2(:,:),   & 
                          EcRPA)

  call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:))
  call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:))

! Compute excitation densities for the T-matrix

  call excitation_density_Tmatrix(ispin,nBas2,nC2,nO2,nR2,nOO,nVV,sERI(:,:,:,:), & 
                                  X1(:,:),Y1(:,:),rho1(:,:,:),             & 
                                  X2(:,:),Y2(:,:),rho2(:,:,:))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  call self_energy_Tmatrix_diag_so(eta,nBas2,nC2,nO2,nV2,nR2,nOO,nVV,seHF(:), & 
                                Omega1(:),rho1(:,:,:),Omega2(:),rho2(:,:,:), & 
                                SigT(:))

! Compute renormalization factor for T-matrix self-energy

  call renormalization_factor_Tmatrix_so(eta,nBas2,nC2,nO2,nV2,nR2,nOO,nVV,seHF(:), & 
                                      Omega1(:),rho1(:,:,:),Omega2(:),rho2(:,:,:), & 
                                      Z(:))

!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------

  eG0T0(:) = seHF(:) + Z(:)*SigT(:)

!----------------------------------------------
! Dump results
!----------------------------------------------

  call print_G0T0(nBas2,nO2,seHF(:),ENuc,ERHF,SigT(:),Z(:),eG0T0(:),EcRPA)

end subroutine soG0T0
