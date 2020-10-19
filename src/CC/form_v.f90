subroutine form_v(nC,nO,nV,nR,X1,X2,X3,X4,t2,v)

! Form quadratic array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: X1(nO-nC,nO-nC,nO-nC,nO-nC)
  double precision,intent(in)   :: X2(nV-nR,nV-nR)
  double precision,intent(in)   :: X3(nO-nC,nO-nC)
  double precision,intent(in)   :: X4(nO-nC,nO-nC,nV-nR,nV-nR)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: v(nO-nC,nO-nC,nV-nR,nV-nR)

  v(:,:,:,:) = 0d0
 
  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR
          do k=1,nO-nC
            do l=1,nO-nC
              v(i,j,a,b) = v(i,j,a,b) + 0.25d0*X1(k,l,i,j)*t2(k,l,a,b)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR
          do c=1,nV-nR
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X2(b,c)*t2(i,j,a,c) + X2(a,c)*t2(i,j,c,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR
          do k=1,nO-nC
            v(i,j,a,b) = v(i,j,a,b) - 0.5d0*(X3(k,j)*t2(i,k,a,b) + X3(k,i)*t2(k,j,a,b))
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR
          do k=1,nO-nC
            do c=1,nV-nR
              v(i,j,a,b) = v(i,j,a,b) + (X4(i,k,a,c)*t2(j,k,b,c) + X4(i,k,b,c)*t2(k,j,a,c))
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_v
