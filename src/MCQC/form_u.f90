subroutine form_u(nO,nV,OOOO,VVVV,OVOV,t2,u)

! Form linear array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: u(nO,nO,nV,nV)

  u(:,:,:,:) = 0d0
 
  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          do c=1,nV
            do d=1,nV
              u(i,j,a,b) = u(i,j,a,b) + 0.5d0*VVVV(a,b,c,d)*t2(i,j,c,d)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do l=1,nO
          do a=1,nV
            do b=1,nV
              u(i,j,a,b) = u(i,j,a,b) + 0.5d0*OOOO(k,l,i,j)*t2(k,l,a,b)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do a=1,nV
          do b=1,nV
            do c=1,nV
              u(i,j,a,b) = u(i,j,a,b) - OVOV(k,b,j,c)*t2(i,k,a,c) &
                                      + OVOV(k,a,j,c)*t2(i,k,b,c) &
                                      - OVOV(k,a,i,c)*t2(j,k,b,c) &
                                      + OVOV(k,b,i,c)*t2(j,k,a,c)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine form_u
