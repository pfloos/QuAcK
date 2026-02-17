subroutine exp_val_1body_rpa(O,V,N,nS,XpY,XmY,Quantity,ExpVal)
! Compute expectation value of a one body quantity in rpa
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)              :: O,V,nS,N
double precision,intent(in)     :: XpY(nS,nS),XmY(nS,nS)
double precision,intent(in)     :: Quantity(N,N)

!Output
double precision,intent(out)    :: ExpVal

!Local
double precision,allocatable    :: XXT(:,:),XYT(:,:),YYT(:,:)
double precision,allocatable    :: X(:,:),Y(:,:)
integer                         :: i,a,jind,b,jb,ia

double precision,external     :: Kronecker_delta

allocate(XXT(nS,nS))
allocate(XYT(nS,nS))
allocate(YYT(nS,nS))
allocate(X(nS,nS))
allocate(Y(nS,nS))

X = 0.5d0*transpose(XpY + XmY)
Y = 0.5d0*transpose(XpY - XmY)
YYT               = matmul(Y,transpose(Y)) 
XYT               = matmul(X,transpose(Y)) 
XXT               = matmul(X,transpose(X))
deallocate(X,Y)

do i = 1, O
  do a = O+1, N
     jind  = i
     b = a
     jb = b - O + (jind - 1) * V 
     ia = a - O +    (i - 1) * V
     
     ! Contributions from fab*dij - fij*dab
     ExpVal = ExpVal + 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))*Quantity(a,b)
     ExpVal = ExpVal - 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))*Quantity(jind,i)
  enddo
enddo
deallocate(XXT,YYT,XYT)
end subroutine
