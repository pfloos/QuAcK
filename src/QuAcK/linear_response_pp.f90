subroutine linear_response_pp(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho,Ec_ppRPA)

! Compute the p-p channel of the linear response: see Scueria et al. JCP 139, 104113 (2013)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  integer                       :: nOO
  integer                       :: nVV
  double precision              :: trace_matrix
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: D(:,:)
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: w1(:)
  double precision,allocatable  :: w2(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

! Output variables

  double precision,intent(out)  :: Ec_ppRPA

! Useful quantities

 nOO = nO*(nO-1)/2
 nVV = nV*(nV-1)/2

! Memory allocation

  allocate(B(nVV,nOO),C(nVV,nVV),D(nOO,nOO),M(nOO+nVV,nOO+nVV),Z(nOO+nVV,nOO+nVV),w(nOO+nVV), &
           w1(nVV),w2(nOO),X1(nVV,nVV),Y1(nVV,nOO),X2(nOO,nVV),Y2(nOO,nOO))

! Build B, C and D matrices for the pp channel

  call linear_response_B_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,B)
  call linear_response_C_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,C)
  call linear_response_D_pp(ispin,dRPA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,D)

!------------------------------------------------------------------------
! Solve the p-p eigenproblem
!------------------------------------------------------------------------
!
!  | C   -B | | X1  X2 |   | w1  0  | | X1  X2 |
!  |        | |        | = |        | |        |
!  | Bt  -D | | Y1  Y2 |   | 0   w2 | | Y1  Y2 |
!

! Diagonal blocks 

  M(    1:nVV    ,    1:nVV)     = + C(1:nVV,1:nVV)
  M(nVV+1:nVV+nOO,nVV+1:nVV+nOO) = - D(1:nOO,1:nOO)

! Off-diagonal blocks

  M(    1:nVV    ,nVV+1:nOO+nVV) = -           B(1:nVV,1:nOO)
  M(nVV+1:nOO+nVV,    1:nVV)     = + transpose(B(1:nVV,1:nOO))

! print*, 'pp-RPA matrix'
! call matout(nOO+nVV,nOO+nVV,M(:,:))

! Diagonalize the p-h matrix

  Z(:,:) = M(:,:)
  call diagonalize_matrix(nOO+nVV,Z(:,:),w(:))

  write(*,*) 'pp-RPA excitation energies'
  call matout(nOO+nVV,1,w(:))
  write(*,*) 

! Split the various quantities in p-p and h-h parts

  w1(:) = w(nOO+1:nOO+nVV)
  w2(:) = w(1:nOO)

  X1(:,:) = Z(nOO+1:nOO+nVV,    1:nVV    )
  Y1(:,:) = Z(nOO+1:nOO+nVV,nVV+1:nOO+nVV)
  X2(:,:) = Z(    1:nOO    ,    1:nVV    ) 
  Y2(:,:) = Z(    1:nOO    ,nVV+1:nOO+nVV)

  if(minval(w1(:)) < 0d0) call print_warning('You may have instabilities in pp-RPA!!')
  if(maxval(w2(:)) > 0d0) call print_warning('You may have instabilities in pp-RPA!!')

! Compute the RPA correlation energy

  Ec_ppRPA = 0.5d0*( sum(w1(:)) - sum(w2(:)) - trace_matrix(nVV,C(:,:)) - trace_matrix(nOO,D(:,:)) )

  print*,'Ec(pp-RPA) = ',Ec_ppRPA
  print*,'Ec(pp-RPA) = ',0.5d0*( sum(abs(w(:))) - trace_matrix(nVV*nOO,M(:,:)))
  print*,'Ec(pp-RPA) = ',+sum(w1(:)) - trace_matrix(nVV,C(:,:))
  print*,'Ec(pp-RPA) = ',-sum(w2(:)) - trace_matrix(nOO,D(:,:))

end subroutine linear_response_pp
