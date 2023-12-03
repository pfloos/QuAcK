subroutine form_u(nC,nO,nV,nR,OOOO,VVVV,OVOV,t2,u)

! Form linear array in CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: OOOO(nO-nC,nO-nC,nO-nC,nO-nC)
  double precision,intent(in)   :: VVVV(nV-nR,nV-nR,nV-nR,nV-nR)
  double precision,intent(in)   :: OVOV(nO-nC,nV-nR,nO-nC,nV-nR)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: u(nO-nC,nO-nC,nV-nR,nV-nR)

  u(:,:,:,:) = 0d0
 
  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR
          do c=1,nV-nR
            do d=1,nV-nR
              u(i,j,a,b) = u(i,j,a,b) + 0.5d0*VVVV(a,b,c,d)*t2(i,j,c,d)
            end do
          end do
        end do
      end do
    end do
  end do

  do i=1,nO-nC
    do j=1,nO-nC
      do k=1,nO-nC
        do l=1,nO-nC
          do a=1,nV-nR
            do b=1,nV-nR
              u(i,j,a,b) = u(i,j,a,b) + 0.5d0*OOOO(k,l,i,j)*t2(k,l,a,b)
            end do
          end do
        end do
      end do
    end do
  end do

  do i=1,nO-nC
    do j=1,nO-nC
      do k=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
            do c=1,nV-nR
              u(i,j,a,b) = u(i,j,a,b) - OVOV(k,b,j,c)*t2(i,k,a,c) &
                                      + OVOV(k,a,j,c)*t2(i,k,b,c) &
                                      - OVOV(k,a,i,c)*t2(j,k,b,c) &
                                      + OVOV(k,b,i,c)*t2(j,k,a,c)
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine 
