subroutine Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY,XmY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb,kc
  double precision              :: delta_spin
  double precision              :: delta_dRPA
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: Ap(:,:)
  double precision,allocatable  :: Bp(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Memory allocation

  allocate(P(nS,nS),Ap(nS,nS),Bp(nS,nS),X(nS,nS),Y(nS,nS))

! Compute P = (X+Y)(X+Y) - 1

  P(:,:) = matmul(transpose(XpY),XpY)

  do ia=1,nS
    P(ia,ia) = P(ia,ia) - 1d0
  enddo

! Compute Aiajb = (ia|bj) and Biajb = (ia|jb)

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

            Ap(ia,jb) = (1d0 + delta_spin)*ERI(i,b,a,j) 
            Bp(ia,jb) = (1d0 + delta_spin)*ERI(i,j,b,a) 

        enddo
      enddo
    enddo
  enddo

! Compute Tr(A x P)

! EcAC = trace_matrix(nS,matmul(Ap,P))

! print*,'EcAC =',EcAC

  X(:,:) = 0.5d0*(XpY(:,:) + XmY(:,:))
  Y(:,:) = 0.5d0*(XpY(:,:) - XmY(:,:))

  EcAC = trace_matrix(nS,matmul(X,matmul(Bp,transpose(Y))) + matmul(Y,matmul(Bp,transpose(X)))) &
       + trace_matrix(nS,matmul(X,matmul(Ap,transpose(X))) + matmul(Y,matmul(Ap,transpose(Y)))) &
       - trace_matrix(nS,Ap)

! print*,'EcAC =',EcAC

end subroutine Ec_AC

