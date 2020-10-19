subroutine form_X(nC,nO,nV,nR,OOVV,t2,X1,X2,X3,X4)

! Form intermediate arrays X's in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: OOVV(nO-nC,nO-nC,nV-nR,nV-nR)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: X1(nO-nC,nO-nC,nO-nC,nO-nC) 
  double precision,intent(out)  :: X2(nV-nR,nV-nR) 
  double precision,intent(out)  :: X3(nO-nC,nO-nC) 
  double precision,intent(out)  :: X4(nO-nC,nO-nC,nV-nR,nV-nR)

! Initialization

  X1(:,:,:,:) = 0d0
  X2(:,:)     = 0d0
  X3(:,:)     = 0d0
  X4(:,:,:,:) = 0d0

! Build X1

  do k=1,nO-nC
    do l=1,nO-nC
      do i=1,nO-nC
        do j=1,nO-nC
          do c=1,nV-nR
            do d=1,nV-nR
              X1(k,l,i,j) = X1(k,l,i,j) + OOVV(k,l,c,d)*t2(i,j,c,d)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X2

  do b=1,nV-nR
    do c=1,nV-nR
      do k=1,nO-nC
        do l=1,nO-nC
          do d=1,nV-nR
            X2(b,c) = X2(b,c) + OOVV(k,l,c,d)*t2(k,l,b,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X3

  do k=1,nO-nC
    do j=1,nO-nC
      do l=1,nO-nC
        do c=1,nV-nR
          do d=1,nV-nR
            X3(k,j) = X3(k,j) + OOVV(k,l,c,d)*t2(j,l,c,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X4

  do i=1,nO-nC
    do l=1,nO-nC
      do a=1,nV-nR
        do d=1,nV-nR
          do k=1,nO-nC
            do c=1,nV-nR
              X4(i,l,a,d) = X4(i,l,a,d) + OOVV(k,l,c,d)*t2(i,k,a,c) 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_X
