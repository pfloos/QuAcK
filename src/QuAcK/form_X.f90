subroutine form_X(nO,nV,OOVV,t2,X1,X2,X3,X4)

! Form intermediate arrays X's in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: X1(nO,nO,nO,nO) 
  double precision,intent(out)  :: X2(nV,nV) 
  double precision,intent(out)  :: X3(nO,nO) 
  double precision,intent(out)  :: X4(nO,nO,nV,nV)

! Initialization

  X1(:,:,:,:) = 0d0
  X2(:,:)     = 0d0
  X3(:,:)     = 0d0
  X4(:,:,:,:) = 0d0

! Build X1

  do k=1,nO
    do l=1,nO
      do i=1,nO
        do j=1,nO
          do c=1,nV
            do d=1,nV
              X1(k,l,i,j) = X1(k,l,i,j) + OOVV(k,l,c,d)*t2(i,j,c,d)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X2

  do b=1,nV
    do c=1,nV
      do k=1,nO
        do l=1,nO
          do d=1,nV
            X2(b,c) = X2(b,c) + OOVV(k,l,c,d)*t2(k,l,b,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X3

  do k=1,nO
    do j=1,nO
      do l=1,nO
        do c=1,nV
          do d=1,nV
            X3(k,j) = X3(k,j) + OOVV(k,l,c,d)*t2(j,l,c,d)
          enddo
        enddo
      enddo
    enddo
  enddo

! Build X4

  do i=1,nO
    do l=1,nO
      do a=1,nV
        do d=1,nV
          do k=1,nO
            do c=1,nV
              X4(i,l,a,d) = X4(i,l,a,d) + OOVV(k,l,c,d)*t2(i,k,a,c) 
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_X
