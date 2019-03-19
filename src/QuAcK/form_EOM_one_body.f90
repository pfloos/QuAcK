subroutine form_EOM_one_body(nO,nV,foo,fov,fvv,t1,t2,OOOV,OOVV,OVVV,cFoo,cFov,cFvv)

! Form one-body part of the EOM matrix

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: foo(nO,nO)
  double precision,intent(in)   :: fov(nO,nV)
  double precision,intent(in)   :: fvv(nV,nV)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)

! Local variables

  integer                       :: i,j,k,l,m,n
  integer                       :: a,b,c,d,e,f
  double precision,allocatable  :: tau(:,:,:,:)

! Output variables

  double precision,intent(out)  :: cFoo(nO,nO)
  double precision,intent(out)  :: cFov(nO,nV)
  double precision,intent(out)  :: cFvv(nV,nV)

  allocate(tau(nO,nO,nV,nV))

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b) - t1(i,b)*t1(j,a)
        end do
      end do
    end do
  end do

! OV block

  cFov(:,:) = 0d0

! VV block

  cFvv(:,:) = fvv(:,:)

  do a=1,nV
    do b=1,nV

      do m=1,nO
        cFvv(a,b) = cFvv(a,b) - fov(m,b)*t1(m,a)
      end do

      do m=1,nO
        do f=1,nV
          cFvv(a,b) = cFvv(a,b) + t1(m,f)*OVVV(m,a,f,b)
        end do
      end do

      do m=1,nO
        do n=1,nO
          do e=1,nV
            cFvv(a,b) = cFvv(a,b) - 0.5d0*tau(m,n,a,e)*OOVV(m,n,b,e)
          end do
        end do
      end do

    enddo
  enddo

! OO block

  cFoo(:,:) = foo(:,:)

  do j=1,nO
    do i=1,nO

      do e=1,nV
        cFoo(j,i) = cFoo(j,i) + t1(i,e)*fov(j,e)
      end do

      do m=1,nO
        do e=1,nV
          cFoo(j,i) = cFoo(j,i) + t1(m,e)*OVVV(j,m,i,e)
        end do
      end do

      do m=1,nO
        do e=1,nV
          do f=1,nV
            cFoo(j,i) = cFoo(j,i) + 0.5d0*tau(i,m,e,f)*OOVV(j,m,e,f)
          end do
        end do
      end do

    end do
  end do

end subroutine form_EOM_one_body
