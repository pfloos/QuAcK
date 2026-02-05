subroutine complex_exp_val_1body_rpa(O,V,N,nS,XpY,XmY,Quantity,ExpVal)
! Compute expectation of a one body quantity in RPA
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)            :: O,V,nS,N
complex*16,intent(in)         :: XpY(nS,nS),XmY(nS,nS)
complex*16,intent(in)         :: Quantity(N,N)

!Output
complex*16,intent(out)        :: ExpVal

!Local
complex*16,allocatable        :: XXT(:,:),XYT(:,:),YYT(:,:)
complex*16,allocatable        :: X(:,:),Y(:,:)
integer                       :: i,a,jind,b,jb,ia

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
