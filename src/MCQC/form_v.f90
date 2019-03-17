subroutine form_v(nO,nV,X1,X2,X3,X4,t2,v)

! Form quadratic array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: X1(nO,nO,nO,nO)
  double precision,intent(in)   :: X2(nV,nV)
  double precision,intent(in)   :: X3(nO,nO)
  double precision,intent(in)   :: X4(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: v(nO,nO,nV,nV)

  v(:,:,:,:) = 0d0
 
  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          do k=1,nO
            do l=1,nO
              v(i,j,a,b) = v(i,j,a,b) + 0.25d0*X1(k,l,i,j)*t2(k,l,a,b)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          do c=1,nV
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X2(b,c)*t2(i,j,a,c) + X2(a,c)*t2(i,j,c,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          do k=1,nO
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X3(k,j)*t2(i,k,a,b) + X3(k,i)*t2(k,j,a,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          do k=1,nO
            do c=1,nV
              v(i,j,a,b) = v(i,j,a,b) + (X4(i,k,a,c)*t2(j,k,b,c) + X4(i,k,b,c)*t2(k,j,a,c))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_v
