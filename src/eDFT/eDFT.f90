program eDFT

! exchange-correlation density-functional theory calculations

  implicit none
  include 'parameters.h'

  integer                       :: nNuc,nBas
  integer                       :: nEl(nspin),nC(nspin),nO(nspin),nV(nspin),nR(nspin)
  double precision              :: ENuc,Ew,EwGIC

  double precision,allocatable  :: ZNuc(:),rNuc(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)
  integer,allocatable           :: max_ang_mom(:)
  double precision,allocatable  :: min_exponent(:,:)
  double precision,allocatable  :: max_exponent(:)

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: ERI(:,:,:,:)
  double precision,allocatable  :: c(:,:)

  character(len=7)              :: method
  integer                       :: x_rung,c_rung
  character(len=12)             :: x_DFA ,c_DFA
  logical                       :: LDA_centered = .true.

  integer                       :: SGn
  double precision              :: radial_precision
  integer                       :: nRad
  integer                       :: nAng
  integer                       :: nGrid
  double precision,allocatable  :: root(:,:)
  double precision,allocatable  :: weight(:)

  double precision,allocatable  :: AO(:,:)
  double precision,allocatable  :: dAO(:,:,:)

  double precision              :: start_KS,end_KS,t_KS
  double precision              :: start_int,end_int,t_int

  integer                       :: nEns
  double precision,allocatable  :: wEns(:)

  integer                       :: maxSCF,max_diis
  double precision              :: thresh
  logical                       :: DIIS
  integer                       :: guess_type
  integer                       :: ortho_type

! Hello World

  write(*,*)
  write(*,*) '******************************************'
  write(*,*) '* eDFT: density-functional for ensembles *'
  write(*,*) '******************************************'
  write(*,*)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electroes of the system
! nC   = number of core orbitals
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nR   = number of Rydberg orbitals 
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nNuc,nEl(:),nO(:),nC(:),nR(:))
  allocate(ZNuc(nNuc),rNuc(nNuc,ncart))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

  allocate(CenterShell(maxShell,ncart),TotAngMomShell(maxShell),KShell(maxShell),DShell(maxShell,maxK), & 
           ExpShell(maxShell,maxK),max_ang_mom(nNuc),min_exponent(nNuc,maxL+1),max_exponent(nNuc))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nNuc,rNuc,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell, & 
                  max_ang_mom,min_exponent,max_exponent)

!------------------------------------------------------------------------
! DFT options
!------------------------------------------------------------------------

! Allocate ensemble weights

  allocate(wEns(maxEns))
  call read_options(method,x_rung,x_DFA,c_rung,c_DFA,SGn, & 
                    nEns,wEns,maxSCF,thresh,DIIS,max_diis,guess_type,ortho_type)

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas), &
           X(nBas,nBas),ERI(nBas,nBas,nBas,nBas),c(nBas,nBas))

 ! Read integrals

  call cpu_time(start_int)

  call system('./GoQCaml')
  call read_integrals(nEl(:),nBas,S,T,V,Hc,ERI)

  call cpu_time(end_int)

  t_int = end_int - start_int
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for reading integrals = ',t_int,' seconds'
  write(*,*)

! Orthogonalization X = S^(-1/2)

  call orthogonalization_matrix(ortho_type,nBas,S,X)

!------------------------------------------------------------------------
! Construct quadrature grid
!------------------------------------------------------------------------
  call read_grid(SGn,radial_precision,nRad,nAng,nGrid)
! nGrid = nRad*nAng

  call allocate_grid(nNuc,ZNuc,max_ang_mom,min_exponent,max_exponent,radial_precision,nAng,nGrid)

  allocate(root(ncart,nGrid),weight(nGrid))

! call quadrature_grid(nRad,nAng,nGrid,root,weight)

  call build_grid(nNuc,ZNuc,rNuc,max_ang_mom,min_exponent,max_exponent, & 
                  radial_precision,nRad,nAng,nGrid,weight,root)

!------------------------------------------------------------------------
! Calculate AO values at grid points
!------------------------------------------------------------------------

  allocate(AO(nBas,nGrid),dAO(ncart,nBas,nGrid))
  call AO_values_grid(nBas,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,nGrid,root,AO,dAO)

  LDA_centered = .true.

!------------------------------------------------------------------------
! Compute GOK-RKS energy
!------------------------------------------------------------------------

  if(method == 'GOK-RKS') then

    call cpu_time(start_KS)
    call GOK_RKS(.false.,x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns(:),nGrid,weight(:), &
                 maxSCF,thresh,max_diis,guess_type,nBas,AO(:,:),dAO(:,:,:),nO(1),nV(1),       &
                 S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,Ew,EwGIC,c(:,:))
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GOK-RKS = ',t_KS,' seconds'
    write(*,*)

 end if

!------------------------------------------------------------------------
! Compute RKS energy
!------------------------------------------------------------------------

  if(method == 'LIM-RKS') then

    call cpu_time(start_KS)
    call LIM_RKS(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns(:),nGrid,weight(:),   &
                 maxSCF,thresh,max_diis,guess_type,nBas,AO(:,:),dAO(:,:,:),nO(1),nV(1), & 
                 S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,c(:,:))
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for LIM-RKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute GOK-UKS energy (BROKEN)
!------------------------------------------------------------------------

  if(method == 'GOK-UKS') then

    call cpu_time(start_KS)
    call GOK_UKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),maxSCF,thresh,max_diis,guess_type, & 
                   nBas,AO(:,:),dAO(:,:,:),nO(:),nV(:),S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,Ew)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! End of eDFT
!------------------------------------------------------------------------
end program eDFT
