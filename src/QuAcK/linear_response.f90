subroutine linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho,EcRPA,Omega,XpY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA,TDA,BSE
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas),rho(nBas,nBas,nS)
  
! Local variables

  double precision              :: trace_matrix
  double precision,allocatable  :: A(:,:),B(:,:),ApB(:,:),AmB(:,:),AmBSq(:,:),Z(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Omega(nS),XpY(nS,nS)


! Memory allocation

  allocate(A(nS,nS),B(nS,nS),ApB(nS,nS),AmB(nS,nS),AmBSq(nS,nS),Z(nS,nS))

! Build A and B matrices 

  call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,e,ERI,A)
  if(BSE) call Bethe_Salpeter_A_matrix(nBas,nC,nO,nV,nR,nS,ERI,Omega,rho,A)

! Tamm-Dancoff approximation

  B = 0d0
  if(.not. TDA) then

    call linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,B)
    if(BSE) call Bethe_Salpeter_B_matrix(nBas,nC,nO,nV,nR,nS,ERI,Omega,rho,B)

  endif

! Build A + B and A - B matrices 

  ApB = A + B
  AmB = A - B

! print*,'A'
! call matout(nS,nS,A)

! print*,'B'
! call matout(nS,nS,B)

! print*,'A+B'
! call matout(nS,nS,ApB)

! print*,'A-B'
! call matout(nS,nS,AmB)

! Diagonalize TD-HF matrix

  call diagonalize_matrix(nS,AmB,Omega)

  if(minval(Omega) < 0d0) &
    call print_warning('You may have instabilities in linear response!!')

  call ADAt(nS,AmB,sqrt(Omega),AmBSq)
  Z = matmul(AmBSq,matmul(ApB,AmBSq))

! print*,'Z'
! call matout(nS,nS,Z)

  call diagonalize_matrix(nS,Z,Omega)

  if(minval(Omega) < 0d0) & 
    call print_warning('You may have instabilities in linear response!!')
 
  Omega = sqrt(Omega)
  XpY = matmul(transpose(Z),AmBSq)
  call DA(nS,1d0/sqrt(abs(Omega)),XpY)

! print*,'X+Y'
! call matout(nS,nS,XpY)

! print*,'RPA excitations'
! call matout(nS,1,Omega)

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Omega) - trace_matrix(nS,A))

! print*,'EcRPA = ',EcRPA

end subroutine linear_response
