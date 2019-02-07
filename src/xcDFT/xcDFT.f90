program xcDFT

! exchange-correlation density-functional theory calculations

  include 'parameters.h'

  integer                       :: nAt,nBas,nEl,nO,nV
  double precision              :: ENuc,EKS

  double precision,allocatable  :: ZNuc(:),rAt(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)

  double precision,allocatable  :: S(:,:),T(:,:),V(:,:),Hc(:,:),X(:,:)
  double precision,allocatable  :: ERI(:,:,:,:)

  integer                       :: rung
  integer                       :: SGn
  integer                       :: nRad,nAng,nGrid
  double precision,allocatable  :: root(:,:)
  double precision,allocatable  :: weight(:)
  double precision,allocatable  :: AO(:,:)
  double precision,allocatable  :: dAO(:,:,:)

  double precision              :: start_KS,end_KS,t_KS

! Hello World

  write(*,*)
  write(*,*) '********************************'
  write(*,*) '* TCCM winter school 2008: DFT *'
  write(*,*) '********************************'
  write(*,*)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nAt,nEl,nO)
  allocate(ZNuc(nAt),rAt(nAt,3))

! Read geometry

  call read_geometry(nAt,ZNuc,rAt,ENuc)

  allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
           DShell(maxShell,maxK),ExpShell(maxShell,maxK))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), &
           ERI(nBas,nBas,nBas,nBas))

!   Read integrals

  call read_integrals(nBas,S,T,V,Hc,ERI)  

! Orthogonalization X = S^(-1/2)

  call orthogonalization_matrix(nBas,S,X)

!------------------------------------------------------------------------
! DFT options
!------------------------------------------------------------------------

  call read_options(rung,SGn)

!------------------------------------------------------------------------
! Construct quadrature grid
!------------------------------------------------------------------------
  call read_grid(SGn,nRad,nAng,nGrid)

  allocate(root(3,nGrid),weight(nGrid))
  call quadrature_grid(nRad,nAng,nGrid,root,weight)

!------------------------------------------------------------------------
! Calculate AO values at grid points
!------------------------------------------------------------------------

  allocate(AO(nBas,nGrid),dAO(3,nBas,nGrid))
  call AO_values_grid(nBas,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      nGrid,root,AO,dAO)

!------------------------------------------------------------------------
! Compute KS energy
!------------------------------------------------------------------------

    call cpu_time(start_KS)
    call RKS(rung,nGrid,weight,nBas,AO,dAO,nO,S,T,V,Hc,ERI,X,ENuc,EKS)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for KS = ',t_KS,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! End of xcDFT
!------------------------------------------------------------------------
end program xcDFT
