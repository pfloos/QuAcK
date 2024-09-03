subroutine pCCD_z_residual(O,V,N,OOVV,OVOV,OVVO,OOOO,VVVV,delta_OV,t2,z2,r2)
      
! Compute the residual for the z amplitudes at the pCCD level

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
  double precision,intent(in)   :: z2(O,V)
      
! Local variables

  integer                       :: i,j,a,b

  double precision,allocatable  :: yO(:,:)
  double precision,allocatable  :: yV(:,:)

! Output variables

  double precision,intent(out)  :: r2(O,V)

! Allocate memory

  allocate(yO(O,O),yV(V,V))

! Form intermediate array

  yO(:,:) = matmul(OOVV,transpose(t2))
  yV(:,:) = matmul(transpose(OOVV),t2)

! Compute residual

  r2(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*z2(:,:) &
          - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - 2d0*OOVV(:,:)*t2(:,:))*z2(:,:)

  do i=1,O
    do a=1,V

      do j=1,O
        r2(i,a) = r2(i,a) - 2d0*OOVV(j,a)*t2(j,a)*z2(i,a) - 2d0*OOVV(i,a)*z2(j,a)*t2(j,a) &
                          + OOOO(i,j)*z2(j,a) + yO(i,j)*z2(j,a)
      end do

      do b=1,V
        r2(i,a) = r2(i,a) - 2d0*OOVV(i,b)*t2(i,b)*z2(i,a) - 2d0*OOVV(i,a)*z2(i,b)*t2(i,b) &
                          + VVVV(b,a)*z2(i,b) + yV(a,b)*z2(i,b)
      end do

    end do
  end do

end
