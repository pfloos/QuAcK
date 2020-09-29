 subroutine MCMP2(doDrift,nBas,nC,nO,nV,c,e,EcMP2,               &
                  nMC,nEq,nWalk,dt,nPrint,                                  &
                  nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                  Norm, &
                  EcMCMP2,Err_EcMCMP2,Var_EcMCMP2)

! Perform Monte Carlo MP2 calculation

  implicit none

  include 'parameters.h'
  include 'quadrature.h'

! Input variables 

  logical,intent(in)            :: doDrift
  integer,intent(in)            :: nBas,nC,nO,nV,nMC,nEq,nWalk,nPrint
  double precision,intent(inout):: dt
  double precision,intent(in)   :: EcMP2(3)
  double precision,intent(in)   :: c(nBas,nBas),e(nBas)

  integer,intent(in)            :: nShell
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: CenterShell(maxShell,3),DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  logical                       :: AcPh,EqPh,Accept,dump
  double precision              :: start_Eq,end_Eq,t_Eq,start_Ac,end_Ac,t_Ac
  integer                       :: nWP
  double precision              :: Norm,NormSq,nData,tau
  double precision,allocatable  :: chi1(:,:,:),chi2(:,:,:),eta(:)
   
  double precision,allocatable  :: cO(:,:),cV(:,:),eO(:),eV(:),P(:,:),eO_Quad(:,:),eV_Quad(:,:)
  double precision,allocatable  :: r(:,:,:), r12(:), gAO(:,:,:),  g(:,:), w(:)
  double precision,allocatable  :: rp(:,:,:),r12p(:),gAOp(:,:,:), gp(:,:),wp(:)
  double precision,allocatable  :: o1MO(:,:),o2MO(:,:),v1MO(:,:),v2MO(:,:)
  double precision,allocatable  :: o11(:,:),o12(:,:),o21(:,:),o22(:,:)
  double precision,allocatable  :: v11(:,:),v12(:,:),v21(:,:),v22(:,:)
  double precision,allocatable  :: fd_Quad(:,:),fx_Quad(:,:),fd(:),fx(:),fdx(:)

  double precision,allocatable  :: dgAO(:,:,:,:),dg(:,:,:),dgAOp(:,:,:,:),dgp(:,:,:)
  double precision,allocatable  :: F(:,:,:),Fp(:,:,:),T(:),Tp(:)

  double precision              :: acceptance,D
  double precision              :: eloc_MP2(3),mean_MP2(3),variance_MP2(3)

  integer                       :: iW,kW,lW,klW,iMC,q

! Output variables

  double precision,intent(out)  :: EcMCMP2(3),Err_EcMCMP2(3),Var_EcMCMP2(3)

! Number of distinct walker pairs

  nWP = nWalk*(nWalk-1)/2

! Diffusion coefficient
 
  D = 0.5d0

! Do diffusion-drift moves?
  
  if(doDrift) then

    write(*,*)
    write(*,*) '*** Diffusion-drift algorithm ***'
    write(*,*)

  else

    write(*,*)
    write(*,*) '*** Diffusion-only  algorithm ***'
    write(*,*)

  endif

! Print results

  dump = .true.
  if(dump) open(unit=13,file='results/data')

!------------------------------------------------------------------------
! Memory allocation
!------------------------------------------------------------------------
  allocate(cO(nBas,nO),cV(nBas,nV),eO(nO),eV(nV),                                   &
           eO_Quad(nQuad,nO),eV_Quad(nQuad,nV),                                     & 
           P(nBas,nBas),r(nWalk,2,3),rp(nWalk,2,3),                                 &
           chi1(nWalk,2,3),chi2(nWalk,2,3),eta(nWalk),                              &
           r12(nWalk),r12p(nWalk),w(nWalk),wp(nWalk),                               &
           g(nWalk,2),gp(nWalk,2),gAO(nWalk,2,nBas),gAOp(nWalk,2,nBas),             &
           dg(nWalk,2,3),dgp(nWalk,2,3),dgAO(nWalk,2,3,nBas),dgAOp(nWalk,2,3,nBas), &
           o1MO(nWalk,nO),v1MO(nWalk,nV),o2MO(nWalk,nO),v2MO(nWalk,nV),             &
           o11(nQuad,nWP),v11(nQuad,nWP),o12(nQuad,nWP),v12(nQuad,nWP),             & 
           o21(nQuad,nWP),v21(nQuad,nWP),o22(nQuad,nWP),v22(nQuad,nWP),             & 
           fd_Quad(nQuad,nWP),fd(nWP),fx_Quad(nQuad,nWP),fx(nWP),fdx(nWP),          &
           T(nWalk),Tp(nWalk),F(nWalk,2,3),Fp(nWalk,2,3))

! Split MOs into occupied and virtual sets

  eO(1:nO) = e(nC+1:nC+nO)
  eV(1:nV) = e(nC+nO+1:nBas)

  do q=1,nQuad
    tau = 1d0/rQuad(q)
    eO_Quad(q,1:nO) = exp(+eO(1:nO)*(tau-1d0))*sqrt(tau)
    eV_Quad(q,1:nV) = exp(-eV(1:nV)*(tau-1d0))*sqrt(tau)
  enddo

  cO(1:nBas,1:nO) = c(1:nBas,nC+1:nC+nO)
  cV(1:nBas,1:nV) = c(1:nBas,nC+nO+1:nBas)

! Compute norm of the trial wave function

  call norm_trial(nBas,nO,cO,P,Norm,NormSq)

!------------------------------------------------------------------------
! Initialize MC-MP2 calculation
!------------------------------------------------------------------------

! Initialize electron coordinates

  call random_number(r)
  r = 2d0*r - 1d0

! Compute initial interelectronic distances

  call rij(nWalk,r,r12)

! Compute initial AO values and their derivatives (if required)

  call AO_values(doDrift,nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,r,gAO,dgAO)

! Compute initial weight function

  call density(doDrift,nBas,nWalk,P,gAO,dgAO,g,dg)

! Compute initial weights

  w(1:nWalk) = g(1:nWalk,1)*g(1:nWalk,2)/r12(1:nWalk)

! Compute initial quantum force

  if(doDrift) call drift(nWalk,r,r12,g,dg,F)

! Equilibration or Accumulation?

  AcPh = .false.
  EqPh = .true.

! Initialization

  nData        = 0d0
  acceptance   = 0d0

  mean_MP2     = 0d0
  variance_MP2 = 0d0
 
  T  = 1d0
  Tp = 1d0
 
!------------------------------------------------------------------------
! Start main Monte Carlo loop
!------------------------------------------------------------------------
  call cpu_time(start_Eq)

  do iMC=1,nEq+nMC

!   Timings 

    if(iMC == nEq + 1) then
      AcPh = .true.
      EqPh = .false.
      write(*,*) 'Time step value at the end of equilibration: dt = ',dt
      write(*,*) 
      call cpu_time(end_Eq)
      t_Eq = end_Eq - start_Eq
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for equilibration = ',t_Eq,' seconds'
      write(*,*)
      call cpu_time(start_Ac)
    endif

!   Optimize time step to reach 50% acceptance

  if(EqPh .and. mod(iMC,100) == 0) call optimize_timestep(nWalk,iMC,acceptance,dt)

!   Move electrons

    call random_number(chi1)
    call random_number(chi2)

!   Diffusion

    rp(:,:,:) = r(:,:,:) + sqrt(2d0*D*dt)*sqrt(-2d0*log(chi1(:,:,:)))*cos(2d0*pi*chi2(:,:,:))

!   Drift 

    if(doDrift) rp(:,:,:) = rp(:,:,:) + D*dt*F(:,:,:)

!   Compute new interelectronic distances

    call rij(nWalk,rp,r12p)

!   Compute new AO values and their derivatives (if required)

    call AO_values(doDrift,nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,rp,gAOp,dgAOp)

    call Density(doDrift,nBas,nWalk,P,gAOp,dgAOp,gp,dgp)

!   Compute new weights

    wp(1:nWalk) = gp(1:nWalk,1)*gp(1:nWalk,2)/r12p(1:nWalk)

!   Compute new quantum force and transition probability
    
    if(doDrift) then

      call Drift(nWalk,rp,r12p,gp,dgp,Fp)
      call transition_probability(nWalk,dt,D,r,rp,F,Fp,T,Tp)

    endif

!   Move for walkers

    call random_number(eta) 

    do iW=1,nWalk

      Accept = (wp(iW)*Tp(iW))/(w(iW)*T(iW)) > eta(iW)

      if(Accept) then

        acceptance = acceptance + 1d0

        r(iW,1:2,1:3) = rp(iW,1:2,1:3)
        gAO(iW,1:2,1:nBas) = gAOp(iW,1:2,1:nBas)
        r12(iW) = r12p(iW)
        w(iW) = wp(iW)

       if(doDrift) F(iW,1:2,1:3) = Fp(iW,1:2,1:3)

      endif

    enddo

!   Accumulation phase 

    if(AcPh) then

      nData = nData + 1d0

!     Calculate Green functions

      call Green_function(nBas,nO,nV,nWalk,nWP,cO,cV,eO_Quad,eV_Quad,gAO, &
                          o1MO,o2MO,v1MO,v2MO,o11,o12,o21,o22,v11,v12,v21,v22)

!     Compute local energy

      fd_Quad = o11*o22*v11*v22 + o12*o21*v12*v21
      fx_Quad = o11*o22*v12*v21 + o12*o21*v11*v22

      fd = matmul(wQuad,fd_Quad)
      fx = matmul(wQuad,fx_Quad)
      
      eloc_MP2 = 0d0
      klW = 0
      do kW=1,nWalk-1
        do lW=kW+1,nWalk
          klW = klW + 1
          eloc_MP2(2) = eloc_MP2(2) + fd(klW)/(r12(kW)*r12(lW)*w(kW)*w(lW))
          eloc_MP2(3) = eloc_MP2(3) + fx(klW)/(r12(kW)*r12(lW)*w(kW)*w(lW))
        enddo
      enddo

      eloc_MP2(2) = -2d0*eloc_MP2(2)/dble(2*nWP)
      eloc_MP2(3) =      eloc_MP2(3)/dble(2*nWP)

      fdx = -2d0*fd + fx
      eloc_MP2(1) = eloc_MP2(2) + eloc_MP2(3)

!     Accumulate results     

      mean_MP2     = mean_MP2     + eloc_MP2 
      variance_MP2 = variance_MP2 + eloc_MP2*eloc_MP2 

!     Print results

      if(mod(iMC,nPrint) == 0) then 

        call compute_error(nData,mean_MP2,variance_MP2,Err_EcMCMP2)
        EcMCMP2     = mean_MP2/nData
        Var_EcMCMP2 = variance_MP2/nData
        EcMCMP2     = Norm*EcMCMP2
        Var_EcMCMP2 = Norm*Var_EcMCMP2
        Err_EcMCMP2 = Norm*Err_EcMCMP2

        write(*,*)
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,I15)')     'Number of data points            ','|',int(nData)
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,10I15)')   'acceptance                       ','|',int(100*acceptance/dble(nWalk*iMC))
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') 'MP2 correlation energy Total     ','|',EcMCMP2(1)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Direct    ','|',EcMCMP2(2)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Exchange  ','|',EcMCMP2(3)
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') 'Statistical error      Total     ','|',Err_EcMCMP2(1)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Direct    ','|',Err_EcMCMP2(2)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Exchange  ','|',Err_EcMCMP2(3)
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') 'Variance               Total     ','|',Var_EcMCMP2(1)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Direct    ','|',Var_EcMCMP2(2)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Exchange  ','|',Var_EcMCMP2(3)
        write(*,*)'-------------------------------------------------------'
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') 'Dev. wrt deterministic Total     ','|',EcMCMP2(1) - EcMP2(1)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Direct    ','|',EcMCMP2(2) - EcMP2(2)
        write(*,'(1X,A36,1X,A1,1X,10F15.8)') '                       Exchange  ','|',EcMCMP2(3) - EcMP2(3)
        write(*,*)'-------------------------------------------------------'
     
        if(dump) write(13,*) int(nData),EcMCMP2(1),Err_EcMCMP2(1)

      endif

    endif

!------------------------------------------------------------------------
! End main Monte Carlo loop
!------------------------------------------------------------------------
  enddo

! Timing 

  call cpu_time(end_Ac)
  t_Ac = end_Ac - start_Ac
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for accumulation = ',t_Ac,' seconds'
  write(*,*)

! Close files

  if(dump)   close(unit=13)

end subroutine MCMP2
