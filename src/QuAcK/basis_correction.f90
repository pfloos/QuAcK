subroutine basis_correction(nBas,nO,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, & 
                            ERI,e,c,P,eG0W0)

! Compute the basis set incompleteness error

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO

  integer,intent(in)            :: nShell
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: CenterShell(maxShell,ncart),DShell(maxShell,maxK),ExpShell(maxShell,maxK)

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)

  double precision,intent(in)   :: eG0W0(nBas)

! Local variables

  integer                       :: SGn
  integer                       :: nRad
  integer                       :: nAng
  integer                       :: nGrid
  double precision,allocatable  :: root(:,:)
  double precision,allocatable  :: weight(:)
  double precision,allocatable  :: AO(:,:)
  double precision,allocatable  :: dAO(:,:,:)
  double precision,allocatable  :: MO(:,:)
  double precision,allocatable  :: dMO(:,:,:)
  double precision,allocatable  :: rho(:)
  double precision,allocatable  :: f(:)
  double precision,allocatable  :: mu(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Basis set incompleteness correction       |'
  write(*,*)'************************************************'
  write(*,*)

! Set quadrature grid

  SGn = 0
  
  call read_grid(SGn,nRad,nAng,nGrid)

! Memory allocation

  allocate(root(ncart,nGrid),weight(nGrid))
  allocate(AO(nBas,nGrid),dAO(ncart,nBas,nGrid),MO(nBas,nGrid),dMO(ncart,nBas,nGrid))
  allocate(rho(nGrid),f(nGrid),mu(nGrid))

  call quadrature_grid(nRad,nAng,nGrid,root,weight)

! Calculate AO values at grid points

  call AO_values_grid(nBas,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,nGrid,root,AO,dAO)

! Calculate MO values at grid points

  call MO_values_grid(nBas,nGrid,c,AO,dAO,MO,dMO)

! Compute one-electron density at grid points 

  call density(nGrid,nBas,P,AO,rho)

! Compute range-sepration function

  call f_grid(nBas,nO,nGrid,MO,ERI,f)
  call mu_grid(nGrid,rho,f,mu)

! Compute energy correction

  call ec_srlda(nGrid,weight,rho,mu)

! Compute orbital corrections

  call fc_srlda(nBas,nGrid,weight,MO,rho,mu,eG0W0)

end subroutine basis_correction
