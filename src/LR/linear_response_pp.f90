subroutine linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOO,nVV, & 
                              e,ERI,Omega1,X1,Y1,Omega2,X2,Y2,EcRPA)

! Compute the p-p channel of the linear response: see Scuseria et al. JCP 139, 104113 (2013)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  
! Local variables

  integer                       :: ab,cd,ij,kl
  integer                       :: p,q,r,s
  double precision              :: trace_matrix
  double precision              :: EcRPA1
  double precision              :: EcRPA2
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
  double precision,intent(out)  :: EcRPA

! Memory allocation

  allocate(B(nVV,nOO),C(nVV,nVV),D(nOO,nOO),M(nOO+nVV,nOO+nVV),Z(nOO+nVV,nOO+nVV),Omega(nOO+nVV))

!-------------------------------------------------!
! Solve the p-p eigenproblem                      !
!-------------------------------------------------!
!                                                 !
!  | C   -B | | X1  X2 |   | w1  0  | | X1  X2 |  !
!  |        | |        | = |        | |        |  !
!  | Bt  -D | | Y1  Y2 |   | 0   w2 | | Y1  Y2 |  !
!                                                 !
!-------------------------------------------------!

! Build B, C and D matrices for the pp channel

  call linear_response_C_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,C)
  call linear_response_D_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,D)

  if(TDA) then

    X1(:,:) = +C(:,:)
    Y1(:,:) = 0d0
    call diagonalize_matrix(nVV,X1,Omega1)
 
    X2(:,:) = 0d0
    Y2(:,:) = -D(:,:)
    call diagonalize_matrix(nOO,Y2,Omega2)

  else

    call linear_response_B_pp(ispin,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,B)

  ! Diagonal blocks 

    M(    1:nVV    ,    1:nVV)     = + C(1:nVV,1:nVV)
    M(nVV+1:nVV+nOO,nVV+1:nVV+nOO) = - D(1:nOO,1:nOO)

  ! Off-diagonal blocks

    M(    1:nVV    ,nVV+1:nOO+nVV) = -           B(1:nVV,1:nOO)
    M(nVV+1:nOO+nVV,    1:nVV)     = + transpose(B(1:nVV,1:nOO))

  ! Diagonalize the p-h matrix

    if(nOO+nVV > 0) call diagonalize_general_matrix(nOO+nVV,M,Omega,Z)

  ! Split the various quantities in p-p and h-h parts

    call sort_ppRPA(nOO,nVV,Omega(:),Z(:,:),Omega1(:),X1(:,:),Y1(:,:),Omega2(:),X2(:,:),Y2(:,:))

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*( sum(Omega1(:)) - sum(Omega2(:)) - trace_matrix(nVV,C(:,:)) - trace_matrix(nOO,D(:,:)) )
  EcRPA1 = +sum(Omega1(:)) - trace_matrix(nVV,C(:,:))
  EcRPA2 = -sum(Omega2(:)) - trace_matrix(nOO,D(:,:))
  if(abs(EcRPA - EcRPA1) > 1d-6 .or. abs(EcRPA - EcRPA2) > 1d-6) & 
    print*,'!!! Issue in pp-RPA linear reponse calculation RPA1 != RPA2 !!!'

end subroutine linear_response_pp
