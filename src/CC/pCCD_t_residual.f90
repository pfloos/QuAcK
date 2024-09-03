subroutine pCCD_t_residual(O,V,N,OOVV,OVOV,OVVO,OOOO,VVVV,delta_OV,t2,r2)
      
! Compute the residual for the t amplitudes at the pCCD level

  implicit none
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  double precision,intent(in)   :: OOVV(O,V)
  double precision,intent(in)   :: OVOV(O,V)
  double precision,intent(in)   :: OVVO(O,V)
  double precision,intent(in)   :: OOOO(O,O)
  double precision,intent(in)   :: VVVV(V,V)
  double precision,intent(in)   :: delta_OV(O,V)
  double precision,intent(in)   :: t2(O,V)
      
! Local variables

  integer                       :: i,j,a,b

  double precision,allocatable  :: yO(:,:)

! Output variables

  double precision,intent(out)  :: r2(O,V)

! Allocate memory

  allocate(yO(O,O))

! Form intermediate array

 yO(:,:) = matmul(t2,transpose(OOVV))

! Compute residual

  r2(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*t2(:,:) &
          - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - OOVV(:,:)*t2(:,:))*t2(:,:)

  do i=1,O
    do a=1,V

      do j=1,O
        r2(i,a) = r2(i,a) - 2d0*OOVV(j,a)*t2(j,a)*t2(i,a) + OOOO(j,i)*t2(j,a) + yO(i,j)*t2(j,a)
      end do

      do b=1,V
        r2(i,a) = r2(i,a) - 2d0*OOVV(i,b)*t2(i,b)*t2(i,a) + VVVV(a,b)*t2(i,b)
      end do

    end do
  end do

end
