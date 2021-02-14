subroutine eDFT(maxSCF,thresh,max_diis,guess_type,mix,nNuc,ZNuc,rNuc,ENuc,nBas,nEl,nC,nO,nV,nR, & 
                nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell, &
                max_ang_mom,min_exponent,max_exponent,S,T,V,Hc,X,ERI,dipole_int,Ew,eKS,cKS,PKS)

! exchange-correlation density-functional theory calculations

! use xc_f90_lib_m

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  logical,intent(in)            :: mix
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nNuc
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEl(nspin)
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)

  integer,intent(in)            :: nShell
  double precision,intent(in)   :: CenterShell(maxShell,ncart)
  integer,intent(in)            :: TotAngMomShell(maxShell)
  integer,intent(in)            :: KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK)
  double precision,intent(in)   :: ExpShell(maxShell,maxK)
  integer,intent(in)            :: max_ang_mom(nNuc)
  double precision,intent(in)   :: min_exponent(nNuc,maxL+1)
  double precision,intent(in)   :: max_exponent(nNuc)


  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  double precision,allocatable  :: c(:,:)

  character(len=8)              :: method
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
  double precision              :: aCC_w1(3)
  double precision              :: aCC_w2(3)

  double precision,allocatable  :: AO(:,:)
  double precision,allocatable  :: dAO(:,:,:)

  double precision              :: start_KS,end_KS,t_KS
  double precision              :: start_int,end_int,t_int

  integer                       :: nEns
  logical                       :: doNcentered
  double precision,allocatable  :: wEns(:)

  double precision,allocatable  :: occnum(:,:,:) 
  integer                       :: Cx_choice

  integer                       :: i,vmajor,vminor,vmicro

! Output variables

  double precision,intent(out)  :: Ew
  double precision,intent(out)  :: eKS(nBas,nspin)
  double precision,intent(out)  :: cKS(nBas,nBas,nspin)
  double precision,intent(out)  :: PKS(nBas,nBas,nspin)
  

! Hello World

  write(*,*)
  write(*,*) '******************************************'
  write(*,*) '* eDFT: density-functional for ensembles *'
  write(*,*) '******************************************'
  write(*,*)

! Libxc version

! call xc_f90_version(vmajor, vminor, vmicro)
! write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

! call xcinfo()

!------------------------------------------------------------------------
! DFT options
!------------------------------------------------------------------------

! Allocate ensemble weights and MO coefficients

  allocate(c(nBas,nspin),wEns(maxEns),occnum(nBas,nspin,maxEns))
  call read_options_dft(nBas,method,x_rung,x_DFA,c_rung,c_DFA,SGn,nEns,wEns,aCC_w1,aCC_w2, & 
                        doNcentered,occnum,Cx_choice)

!------------------------------------------------------------------------
! Construct quadrature grid
!------------------------------------------------------------------------
  call read_grid(SGn,radial_precision,nRad,nAng,nGrid)

  call allocate_grid(nNuc,ZNuc,max_ang_mom,min_exponent,max_exponent,radial_precision,nAng,nGrid)

  allocate(root(ncart,nGrid),weight(nGrid))

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
    call GOK_RKS(.false.,x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight, &
                 maxSCF,thresh,max_diis,guess_type,nBas,AO,dAO,nO(1),nV(1), &
                 S,T,V,Hc,ERI,X,ENuc,Ew,c,occnum,Cx_choice)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for GOK-RKS = ',t_KS,' seconds'
    write(*,*)

 end if

!------------------------------------------------------------------------
! Compute LIM excitation energies
!------------------------------------------------------------------------

  if(method == 'LIM-RKS') then

    call cpu_time(start_KS)
    call LIM_RKS(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,nGrid,weight(:),           &
                 aCC_w1,aCC_w2,maxSCF,thresh,max_diis,guess_type,nBas,AO(:,:),dAO(:,:,:),nO(1),nV(1), & 
                 S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,c(:,:),occnum,Cx_choice,doNcentered)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for LIM-RKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute MOM excitation energies
!------------------------------------------------------------------------

  if(method == 'MOM-RKS') then

    call cpu_time(start_KS)
    call MOM_RKS(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,nGrid,weight(:),           &
                 aCC_w1,aCC_w2,maxSCF,thresh,max_diis,guess_type,nBas,AO(:,:),dAO(:,:,:),nO(1),nV(1), & 
                 S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,c(:,:),occnum,Cx_choice,doNcentered)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MOM-RKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute GOK-UKS energy (BROKEN)
!------------------------------------------------------------------------

  if(method == 'GOK-UKS') then

    call cpu_time(start_KS)
    call GOK_UKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),aCC_w1,aCC_w2,maxSCF,thresh,max_diis,guess_type, & 
                   nBas,AO(:,:),dAO(:,:,:),nO(:),nV(:),S(:,:),T(:,:),V(:,:),Hc(:,:),ERI(:,:,:,:),X(:,:),ENuc,Ew,occnum, &
                   Cx_choice,doNcentered)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute N-centered UKS energy 
!------------------------------------------------------------------------

  if(method == 'eDFT-UKS') then

    call cpu_time(start_KS)
    call eDFT_UKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight(:),maxSCF,thresh,max_diis,guess_type,mix, & 
                  nBas,AO,dAO,S,T,V,Hc,ERI,X,ENuc,occnum,Cx_choice,doNcentered,Ew,eKS,cKS,PKS)
    call cpu_time(end_KS)

    t_KS = end_KS - start_KS
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UKS = ',t_KS,' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! End of eDFT
!------------------------------------------------------------------------
end subroutine eDFT
