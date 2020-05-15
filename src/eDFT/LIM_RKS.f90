subroutine LIM_RKS(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,nGrid,weight, &
                   maxSCF,thresh,max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,c)

! Perform restricted Kohn-Sham calculation for ensembles

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: maxSCF,max_diis,guess_type
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(out)  :: c(nBas,nBas) 

! Local variables

  integer                       :: iEns
  double precision              :: Ew(nEnS)
  double precision              :: wLIM(nEns)
  double precision              :: Om(nEns)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*        Linear-interpolation method           *'
  write(*,*)'*          for excitation energies             *'
  write(*,*)'************************************************'
  write(*,*)

! Initializatio

  Ew(:) = 0d0
  Om(:) = 0d0

!------------------------------------------------------------------------
! Zero-weight calculation
!------------------------------------------------------------------------

  write(*,'(A40)')           '*************************************************'
  write(*,'(A40)')           '            ZERO-WEIGHT CALCULATION              '
  write(*,'(A40)')           '*************************************************'

  wLIM(1) = 1d0
  wLIM(2) = 0d0
  wLIM(3) = 0d0

  do iEns=1,nEns
    write(*,'(A20,I2,A2,F16.10)') ' Weight of state  ',iEns,': ',wLIM(iEns)
  end do
  write(*,'(A40)')           '*************************************************'
  write(*,*)

  call GOK_RKS(.false.,x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wLIM,nGrid,weight,maxSCF,thresh, &
               max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,Ew(1),c)

!------------------------------------------------------------------------
! Equiensemble calculation
!------------------------------------------------------------------------

  write(*,'(A40)')           '*************************************************'
  write(*,'(A40)')           '  TWO-STATE   EQUI-WEIGHT CALCULATION            '
  write(*,'(A40)')           '*************************************************'

  wLIM(1) = 0.5d0
  wLIM(2) = 0.0d0
  wLIM(3) = 0.5d0

  do iEns=1,nEns
    write(*,'(A20,I2,A2,F16.10)') ' Weight of state  ',iEns,': ',wLIM(iEns)
  end do
  write(*,'(A40)')           '*************************************************'
  write(*,*)

  call GOK_RKS(.true.,x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wLIM,nGrid,weight,maxSCF,thresh, &
               max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,Ew(2),c)

!------------------------------------------------------------------------
! Equiensemble calculation
!------------------------------------------------------------------------

  write(*,'(A40)')           '*************************************************'
  write(*,'(A40)')           '  THREE-STATE EQUI-WEIGHT CALCULATION            '
  write(*,'(A40)')           '*************************************************'

  wLIM(1) = 1d0/3d0
  wLIM(2) = 1d0/3d0
  wLIM(3) = 1d0/3d0

  do iEns=1,nEns
    write(*,'(A20,I2,A2,F16.10)') ' Weight of state  ',iEns,': ',wLIM(iEns)
  end do
  write(*,'(A40)')           '*************************************************'
  write(*,*)

! call GOK_RKS(.true.,x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wLIM,nGrid,weight,maxSCF,thresh, &
!              max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,Ew(3),c)


!------------------------------------------------------------------------
! LIM excitation energies
!------------------------------------------------------------------------

  Om(2) = 2d0*(Ew(2) - Ew(1))
! Om(3) = 3d0*(Ew(3) - Ew(2)) + 0.5d0*Om(2)

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' LINEAR INTERPOLATION METHOD EXCITATION ENERGIES '
  write(*,'(A60)')           '-------------------------------------------------'

  write(*,'(A44,F16.10,A3)') '     Ensemble   energy  #1      ',Ew(1),' au'
  write(*,'(A44,F16.10,A3)') '     Ensemble   energy  #2      ',Ew(2),' au'
  write(*,'(A44,F16.10,A3)') '     Ensemble   energy  #3      ',Ew(3),' au'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') '    Excitation energy  1 ->',iEns,': ',Om(iEns),       ' au'
  end do
  write(*,*)
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') '    Excitation energy  1 ->',iEns,': ',Om(iEns)*HaToeV,' eV'
  end do
  write(*,'(A60)')           '-------------------------------------------------'


end subroutine LIM_RKS
