subroutine linear_response_pp(ispin,BSE,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,Omega1,X1,Y1,Omega2,X2,Y2,Ec_ppRPA)

! Compute the p-p channel of the linear response: see Scueria et al. JCP 139, 104113 (2013)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  
! Local variables

  double precision              :: trace_matrix
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: D(:,:)
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: Omega(:)

! Output variables

  double precision,intent(out)  :: Omega1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Omega2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)
  double precision,intent(out)  :: Ec_ppRPA

! Memory allocation

  allocate(B(nVV,nOO),C(nVV,nVV),D(nOO,nOO),M(nOO+nVV,nOO+nVV),Z(nOO+nVV,nOO+nVV),Omega(nOO+nVV))

! Build B, C and D matrices for the pp channel

  call linear_response_B_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,B)
  call linear_response_C_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,C)
  call linear_response_D_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,D)

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
  call diagonalize_matrix(nOO+nVV,Z(:,:),Omega(:))

! write(*,*) 'pp-RPA excitation energies'
! call matout(nOO+nVV,1,Omega(:))
! write(*,*) 

! Split the various quantities in p-p and h-h parts

  Omega1(:) = Omega(nOO+1:nOO+nVV)
  Omega2(:) = Omega(1:nOO)

  X1(:,:) = Z(nOO+1:nOO+nVV,nOO+1:nOO+nVV)
  Y1(:,:) = Z(    1:nOO    ,nOO+1:nOO+nVV)
  X2(:,:) = Z(nOO+1:nOO+nVV,    1:nOO    ) 
  Y2(:,:) = Z(    1:nOO    ,nOO+1:nOO+nVV)

  if(minval(Omega1(:)) < 0d0) call print_warning('You may have instabilities in pp-RPA!!')
  if(maxval(Omega2(:)) > 0d0) call print_warning('You may have instabilities in pp-RPA!!')

! Compute the RPA correlation energy

  Ec_ppRPA = 0.5d0*( sum(Omega1(:)) - sum(Omega2(:)) - trace_matrix(nVV,C(:,:)) - trace_matrix(nOO,D(:,:)) )

  print*,'Ec(pp-RPA) = ',Ec_ppRPA
  print*,'Ec(pp-RPA) = ',+sum(Omega1(:)) - trace_matrix(nVV,C(:,:))
  print*,'Ec(pp-RPA) = ',-sum(Omega2(:)) - trace_matrix(nOO,D(:,:))

end subroutine linear_response_pp
