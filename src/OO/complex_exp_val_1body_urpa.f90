subroutine complex_exp_val_1body_urpa(nO,nV,nBas,nSa,nSb,nSt,nFC,nCVS,XpY,XmY,occupations,virtuals,Quantity,ExpVal)
! Compute expectation of a one body quantity in RPA
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)            :: nO(nspin),nV(nspin),nSa,nSb,nSt,nBas,nFC(nspin),nCVS(nspin)
complex*16,intent(in)         :: XpY(nSt,nSt),XmY(nSt,nSt)
complex*16,intent(in)         :: Quantity(nBas,nBas)
integer,intent(in)            :: occupations(maxval(nO),nspin)
integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)

!Output
complex*16,intent(out)        :: ExpVal

!Local
complex*16,allocatable        :: XXT(:,:),XYT(:,:),YYT(:,:)
complex*16,allocatable        :: X(:,:),Y(:,:)
integer                       :: i,a,j,b,jb,ia

double precision,external     :: Kronecker_delta

allocate(XXT(nSt,nSt))
allocate(XYT(nSt,nSt))
allocate(YYT(nSt,nSt))
allocate(X(nSt,nSt))
allocate(Y(nSt,nSt))

X = 0.5d0*transpose(XpY + XmY)
Y = 0.5d0*transpose(XpY - XmY)
YYT               = matmul(Y,transpose(Y)) 
XYT               = matmul(X,transpose(Y)) 
XXT               = matmul(X,transpose(X))
deallocate(X,Y)

! aaaa block
ia = 0
do i=1,nO(1) - nFC(1)
  do a=nCVS(1)+1,nBas - nO(1)
    ia = ia + 1
    jb = 0
    do j=1,nO(1) - nFC(1)
      do b=nCVS(1)+1,nBas - nO(1)
        jb = jb + 1
    
        if(i==j) then
          ExpVal = ExpVal + 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))*Quantity(virtuals(a,1),virtuals(b,1))
        endif
        if(a==b) then
          ExpVal = ExpVal - 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))*Quantity(occupations(j,1),occupations(i,1))
        endif

      enddo
    enddo
  enddo
enddo

! bbbb block
ia = 0
do i=1,nO(1) - nFC(1)
  do a=nCVS(1)+1,nBas - nO(1)
    ia = ia + 1
    jb = 0
    do j=1,nO(1) - nFC(1)
      do b=nCVS(1)+1,nBas - nO(1)
        jb = jb + 1
    
        if(i==j) then
          ExpVal = ExpVal + 0.5d0*(XXT(nSa + jb,nSa + ia) + YYT(nSa + jb,nSa + ia) - Kronecker_delta(jb,ia))*Quantity(virtuals(a,1),virtuals(b,1))
        endif
        if(a==b) then
          ExpVal = ExpVal - 0.5d0*(XXT(nSa + jb,nSa + ia) + YYT(nSa + jb,nSa + ia) - Kronecker_delta(jb,ia))*Quantity(occupations(j,2),occupations(i,2))
        endif

      enddo
    enddo
  enddo
enddo
deallocate(XXT,YYT,XYT)
end subroutine
