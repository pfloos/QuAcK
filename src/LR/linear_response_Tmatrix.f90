subroutine linear_response_Tmatrix(ispin,dRPA,TDA,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A_BSE,B_BSE,EcRPA,Omega,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA,TDA
  double precision,intent(in)   :: eta
  integer,intent(in)            :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: A_BSE(nS,nS)
  double precision,intent(in)   :: B_BSE(nS,nS)

! Local variables

  double precision              :: trace_matrix
  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Omega(nS)
  double precision,intent(out)  :: XpY(nS,nS)
  double precision,intent(out)  :: XmY(nS,nS)

! Memory allocation

  allocate(A(nS,nS),B(nS,nS),ApB(nS,nS),AmB(nS,nS),AmBSq(nS,nS),AmBIv(nS,nS),Z(nS,nS))

! Build A and B matrices 

  call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,A)

  print*,'A'
  call matout(nS,nS,A)
  print*,'TA'
  call matout(nS,nS,A_BSE)

  A(:,:) = A(:,:) + A_BSE(:,:)
! if(ispin == 1) A(:,:) = A(:,:) + A_BSE(:,:)
! if(ispin == 2) A(:,:) = A(:,:) - A_BSE(:,:)

! Tamm-Dancoff approximation

  if(TDA) then
 
    B(:,:)   = 0d0
    XpY(:,:) = A(:,:)
    call diagonalize_matrix(nS,XpY,Omega)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    call linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,B)

    print*,'B'
    call matout(nS,nS,B)
    print*,'TB'
    call matout(nS,nS,B_BSE)

    B(:,:) = B(:,:) + B_BSE(:,:)
!   if(ispin == 1) B(:,:) = B(:,:) + B_BSE(:,:)
!   if(ispin == 2) B(:,:) = B(:,:) - B_BSE(:,:)

    ! Build A + B and A - B matrices 

    ApB = A + B
    AmB = A - B

  ! Diagonalize linear response matrix

    call diagonalize_matrix(nS,AmB,Omega)

    if(minval(Omega) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')

!   do ia=1,nS
!     if(Omega(ia) < 0d0) Omega(ia) = 0d0
!   end do

    call ADAt(nS,AmB,1d0*sqrt(Omega),AmBSq)
    call ADAt(nS,AmB,1d0/sqrt(Omega),AmBIv)

    Z = matmul(AmBSq,matmul(ApB,AmBSq))

   call diagonalize_matrix(nS,Z,Omega)

    if(minval(Omega) < 0d0) & 
      call print_warning('You may have instabilities in linear response: negative excitations!!')
 
  ! do ia=1,nS
  !   if(Omega(ia) < 0d0) Omega(ia) = 0d0
  ! end do

    Omega = sqrt(Omega)

    XpY = matmul(transpose(Z),AmBSq)
    call DA(nS,1d0/sqrt(Omega),XpY)

    XmY = matmul(transpose(Z),AmBIv)
    call DA(nS,1d0*sqrt(Omega),XmY)
 
  end if

  ! Compute the RPA correlation energy

    EcRPA = 0.5d0*(sum(Omega) - trace_matrix(nS,A))

end subroutine linear_response_Tmatrix
