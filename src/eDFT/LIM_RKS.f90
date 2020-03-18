subroutine LIM_RKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nGrid,weight,maxSCF,thresh, & 
                   max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,F)

! Perform restricted Kohn-Sham calculation for ensembles

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
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

  double precision,intent(out)  :: F(nBas,nBas) 

! Local variables

  integer                       :: iEns
  double precision              :: EwZW,EwGICZW
  double precision              :: EwEW,EwGICEW
  double precision              :: wLIM(nEns)
  double precision              :: Om(nEns),OmGIC(nEns)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*        Linear-interpolation method           *'
  write(*,*)'*          for excitation energies             *'
  write(*,*)'************************************************'
  write(*,*)

!------------------------------------------------------------------------
! Zero-weight calculation
!------------------------------------------------------------------------

  write(*,'(A40)')           '*************************************************'
  write(*,'(A40)')           '            ZERO-WEIGHT CALCULATION              '
  write(*,'(A40)')           '*************************************************'

  wLIM(1)      = 1d0
  wLIM(2:nEns) = 0d0

  do iEns=1,nEns
    write(*,'(A20,I2,A2,F16.10)') ' Weight of state  ',iEns,': ',wLIM(iEns)
  end do
  write(*,'(A40)')           '*************************************************'
  write(*,*)

  call GOK_RKS(.false.,x_rung,x_DFA,c_rung,c_DFA,nEns,wLIM,nGrid,weight,maxSCF,thresh, &
               max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,EwZW,EwGICZW,F)

!------------------------------------------------------------------------
! Equiensemble calculation
!------------------------------------------------------------------------

  write(*,'(A40)')           '*************************************************'
  write(*,'(A40)')           '            ZERO-WEIGHT CALCULATION              '
  write(*,'(A40)')           '*************************************************'

  wLIM(1:nEns) = 1d0/dble(nEns)

  do iEns=1,nEns
    write(*,'(A20,I2,A2,F16.10)') ' Weight of state  ',iEns,': ',wLIM(iEns)
  end do
  write(*,'(A40)')           '*************************************************'
  write(*,*)

  call GOK_RKS(.true.,x_rung,x_DFA,c_rung,c_DFA,nEns,wLIM,nGrid,weight,maxSCF,thresh, &
               max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,EwEW,EwGICEW,F)

!------------------------------------------------------------------------
! LIM excitation energies
!------------------------------------------------------------------------

  Om(1)    = 0d0
  Om(2)    = 2d0*(EwEW    - EwZW)

  OmGIC(1) = 0d0
  OmGIC(2) = 2d0*(EwGICEW - EwGICZW)

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' LINEAR INTERPOLATION METHOD EXCITATION ENERGIES '
  write(*,'(A60)')           '-------------------------------------------------'

  write(*,'(A44,F16.10,A3)') ' Zero-weight     ensemble energy',EwZW,   ' au'
  write(*,'(A44,F16.10,A3)') ' Zero-weight GIC ensemble energy',EwGICZW,' au'
  write(*,*)
  write(*,'(A44,F16.10,A3)') ' Equi-weight     ensemble energy',EwEW,    ' au'
  write(*,'(A44,F16.10,A3)') ' Equi-weight GIC ensemble energy',EwGICEW, ' au'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') '    Excitation energy  1 ->',iEns,': ',Om(iEns),       ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') '    Excitation energy  1 ->',iEns,': ',Om(iEns)*HaToeV,' eV'
  end do
  write(*,*)
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' GIC Excitation energy  1 ->',iEns,': ',OmGIC(iEns),       ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') ' GIC Excitation energy  1 ->',iEns,': ',OmGIC(iEns)*HaToeV,' eV'
  end do
  write(*,'(A60)')           '-------------------------------------------------'


end subroutine LIM_RKS
