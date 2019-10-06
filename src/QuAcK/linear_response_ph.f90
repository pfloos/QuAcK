subroutine linear_response_ph(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho,Ec_phRPA,Omega,XpY)

! Compute the p-h channel of the linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA,TDA,BSE
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas),rho(nBas,nBas,nS)
  
! Local variables

  double precision,external     :: trace_matrix
  double precision,allocatable  :: A(:,:),B(:,:),M(:,:),w(:)

! Output variables

  double precision,intent(out)  :: Ec_phRPA
  double precision,intent(out)  :: Omega(nS),XpY(nS,nS)


! Memory allocation

  allocate(A(nS,nS),B(nS,nS),M(2*nS,2*nS),w(2*nS))

! Build A and B matrices 

  call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,e,ERI,A)
  if(BSE) call Bethe_Salpeter_A_matrix(nBas,nC,nO,nV,nR,nS,ERI,Omega,rho,A)

! Tamm-Dancoff approximation

  B(:,:) = 0d0
  if(.not. TDA) then

    call linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,B)
    if(BSE) call Bethe_Salpeter_B_matrix(nBas,nC,nO,nV,nR,nS,ERI,Omega,rho,B)

  endif

!------------------------------------------------------------------------
! Solve the p-h eigenproblem
!------------------------------------------------------------------------
!
!  | +A   +B | | X  Y |   | w   0 | | X  Y | 
!  |         | |      | = |       | |      | 
!  | -B   -A | | Y  X |   | 0  -w | | Y  X | 
!

! Diagonal blocks 

  M(1:nS,1:nS)           = +A(1:nS,1:nS)
  M(nS+1:2*nS,nS+1:2*nS) = -A(1:nS,1:nS)


! Off-diagonal blocks

  M(1:nS,nS+1:2*nS) = -B(1:nS,1:nS)
  M(nS+1:2*nS,1:nS) = +B(1:nS,1:nS)

! Diagonalize the p-h matrix

  call diagonalize_matrix(2*nS,M(:,:),w(:))

  Omega(1:nS) = w(nS+1:2*nS)

! Build X+Y

  XpY(1:nS,1:nS) = M(nS+1:2*nS,1:nS) + M(nS+1:2*nS,nS+1:2*nS)

  call DA(nS,1d0/sqrt(Omega),XpY)

! print*,'X+Y'
! call matout(nS,nS,XpY)

! print*,'RPA excitations'
  call matout(2*nS,1,w(:))

! Compute the RPA correlation energy

  Ec_phRPA = 0.5d0*(sum(Omega) - trace_matrix(nS,A))

  print*,'Ec_phRPA = ',Ec_phRPA

end subroutine linear_response_ph
