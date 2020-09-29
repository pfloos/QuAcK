subroutine MinMCMP2(nBas,nEl,nC,nO,nV,c,e,EcMP2,                              &
                    nMC,nEq,nWalk,dt,nPrint,                                  &
                    nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                    TrialType,Norm,cTrial,gradient,hessian)

! Minimize the variance of MC-MP2 

  implicit none
  include 'parameters.h'

! Input variables 

  integer,intent(in)            :: nBas,nEl,nC,nO,nV,nMC,nEq,nWalk,nPrint
  double precision,intent(in)   :: EcMP2(3),dt
  double precision,intent(in)   :: c(nBas,nBas),e(nBas)

  integer,intent(in)            :: nShell
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: CenterShell(maxShell,3),DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  logical                       :: debug,varmin,mincvg
  double precision              :: thresh
  double precision,allocatable  :: max_gradient(:),energy_MCMP2(:),variance_MCMP2(:),error_MCMP2(:),NormTr(:)

  double precision              :: EcMCMP2(3),Err_EcMCMP2(3),Var_EcMCMP2(3)

  integer                       :: it,nIt,i

! Output variables

  integer,intent(in)            :: TrialType
  double precision,intent(inout):: Norm,cTrial(nBas),gradient(nBas),hessian(nBas,nBas)

! Debuging mode

!  debug = .true.
  debug = .false.

! Minimization parameters

  varmin = .true.
  mincvg = .false.
  nIt = 10
  thresh = 1d-5
  allocate(max_gradient(nIt),energy_MCMP2(nIt),variance_MCMP2(nIt),error_MCMP2(nIt),normTr(nIt))

  if(TrialType == 1) then

!   Use HOMO as guess
    cTrial(1:nBas) = c(1:nBas,nEl/2)
!   Normalization factor will be computed later

  endif

!------------------------------------------------------------------------
! Start MC-MP2 variance minimization
!------------------------------------------------------------------------
  it = 0
  do while (it < nIt .and. .not.(mincvg))

    it = it + 1

    write(*,*) '**********************************************************************'
    write(*,*) '               Variance minimization of MC-MP2 energy                 '
    write(*,*) '**********************************************************************'
    write(*,*) ' Iteration n.',it
    write(*,*) '**********************************************************************'

    write(*,*)
    write(*,*) ' Trial wave function coefficients at iteration n.',it 
    call matout(nBas,1,cTrial)
    write(*,*)

    call MCMP2(varmin,nBas,nEl,nC,nO,nV,c,e,EcMP2,                 &
               nMC,nEq,nWalk,dt,nPrint,                                  &
               nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
               TrialType,Norm,cTrial,gradient,hessian,            &
               EcMCMP2,Err_EcMCMP2,Var_EcMCMP2)

!   Newton update of the coefficients

    call Newton(nBas,gradient,hessian,cTrial)

!   Check for convergence
    
    max_gradient(it)   = maxval(abs(gradient))
    energy_MCMP2(it)   = EcMCMP2(1)
    variance_MCMP2(it) = Var_EcMCMP2(1)
    error_MCMP2(it)    = Err_EcMCMP2(1)
    NormTr(it)         = Norm

    write(*,*)
    write(*,*) 'Maximum gradient at iteration n.',it,':',max_gradient(it)
    write(*,*)

    if(max_gradient(it) < thresh) then
      write(*,*) ' Miracle! Variance minimization of MC-MP2 has converged!'
      mincvg = .true.
    endif
 
  enddo

  write(*,*) 
  write(*,*) '********************************'
  write(*,*) 'Summary of variance minimization'
  write(*,*) '********************************'
  write(*,*) 

  write(*,'(A3,A20,A20,A20,A20,A20,A20)') & 
    'It.','Gradient','Ec(MC-MPC2)','Variance','Error','Ec(MC-MP2)-Ec(MP2)','Norm'
  write(*,'(I3,4X,F16.10,4X,F16.10,4X,F16.10,4X,F16.10,4X,F16.10,4X,F16.10)') & 
    (i,max_gradient(i),energy_MCMP2(i),variance_MCMP2(i),error_MCMP2(i),energy_MCMP2(i)-EcMP2(1),NormTr(i),i=1,it)
  write(*,*) 

!------------------------------------------------------------------------
! End MC-MP2 variance minimization
!------------------------------------------------------------------------

end subroutine MinMCMP2
