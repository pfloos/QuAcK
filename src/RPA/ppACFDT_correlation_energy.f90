subroutine ppACFDT_correlation_energy(ispin,nBas,nC,nO,nV,nR,nS,ERI,nOO,nVV,X1,Y1,X2,Y2,EcAC)

! Compute the correlation energy via the adiabatic connection formula for the pp sector

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: D(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC
  
! Memory allocation

  allocate(B(nVV,nOO),C(nVV,nVV),D(nOO,nOO))

! Build pp matrices

  call ppLR_B   (ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,B)
  call ppLR_C_od(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,C)
  call ppLR_D_od(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,D)

! Compute Tr(K x P_lambda)

  EcAC = trace_matrix(nVV,matmul(transpose(X1),matmul(B,Y1)) + matmul(transpose(Y1),matmul(transpose(B),X1))) &
       + trace_matrix(nVV,matmul(transpose(X1),matmul(C,X1)) + matmul(transpose(Y1),matmul(D,Y1)))  &
       + trace_matrix(nOO,matmul(transpose(X2),matmul(B,Y2)) + matmul(transpose(Y2),matmul(transpose(B),X2)))  &
       + trace_matrix(nOO,matmul(transpose(X2),matmul(C,X2)) + matmul(transpose(Y2),matmul(D,Y2)))  &
       - trace_matrix(nVV,C) - trace_matrix(nOO,D)

end subroutine 

