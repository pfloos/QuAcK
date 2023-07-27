subroutine form_EOM_two_body(nO,nV,foo,fov,fvv,t1,t2,cFoo,cFov,cFvv,                      & 
                             OOOO,VOOO,OVOO,OOVO,OOOV,OOVV,VOOV,OVVO,OVVV,VOVV,VVVO,VVVV, & 
                             cWvvoo,cWoooo,cWvvvv,cWvovv,cWooov,cWoovv,cWvoov,cWvvvo,cWovoo)

! Form two-body part of the EOM matrix

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: foo(nO,nO)
  double precision,intent(in)   :: fov(nO,nV)
  double precision,intent(in)   :: fvv(nV,nV)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: VOOO(nV,nO,nO,nO)
  double precision,intent(in)   :: OVOO(nO,nV,nO,nO)
  double precision,intent(in)   :: OOVO(nO,nO,nV,nO)
  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: VOOV(nV,nO,nO,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)
  double precision,intent(in)   :: VOVV(nV,nO,nV,nV)
  double precision,intent(in)   :: VVVO(nV,nV,nV,nO)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)

  double precision,intent(in)   :: cFoo(nO,nO)
  double precision,intent(in)   :: cFov(nO,nV)
  double precision,intent(in)   :: cFvv(nV,nV)

! Local variables

  integer                       :: i,j,k,l,m,n
  integer                       :: a,b,c,d,e,f
  double precision,allocatable  :: tau(:,:,:,:)

! Output variables

  double precision,intent(out)  :: cWvvoo(nV,nV,nO,nO)
  double precision,intent(out)  :: cWoooo(nO,nO,nO,nO)
  double precision,intent(out)  :: cWvvvv(nV,nV,nV,nV)
  double precision,intent(out)  :: cWvovv(nV,nO,nV,nV)
  double precision,intent(out)  :: cWooov(nO,nO,nO,nV)
  double precision,intent(out)  :: cWoovv(nO,nO,nV,nV)
  double precision,intent(out)  :: cWvoov(nV,nO,nO,nV)
  double precision,intent(out)  :: cWvvvo(nV,nV,nV,nO)
  double precision,intent(out)  :: cWovoo(nO,nV,nO,nO)

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

! VVOO block

  cWvvoo(:,:,:,:) = 0d0

! OOOO block

  cWoooo(:,:,:,:) = OOOO(:,:,:,:)

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do l=1,nO

          do e=1,nV
            cWoooo(i,j,k,l) = cWoooo(i,j,k,l) + t1(j,e)*OOOV(i,j,k,e)
            cWoooo(i,j,k,l) = cWoooo(i,j,k,l) - t1(i,e)*OOOV(i,j,l,e)
          end do

          do a=1,nV
            do b=1,nV
              cWoooo(i,j,k,l) = cWoooo(i,j,k,l) + 0.5d0*tau(k,l,a,b)*OOVV(i,j,a,b)
            end do
          end do

        end do
      end do
    end do
  end do

! VVVV block

  cWvvvv(:,:,:,:) = VVVV(:,:,:,:)

  do a=1,nV
    do b=1,nV
      do c=1,nV
        do d=1,nV

          do i=1,nO
            cWvvvv(a,b,c,d) = cWvvvv(a,b,c,d) - t1(i,b)*VOVV(a,i,c,d)
            cWvvvv(a,b,c,d) = cWvvvv(a,b,c,d) + t1(i,a)*VOVV(b,i,c,d)
          end do

          do i=1,nO
            do j=1,nO
              cWvvvv(a,b,c,d) = cWvvvv(a,b,c,d) + tau(i,j,a,b)*OOVV(i,j,c,d)
            end do
          end do

        end do
      end do
    end do
  end do

! VOVV block

  cWvovv(:,:,:,:) = VOVV(:,:,:,:)

  do a=1,nV
    do i=1,nO
      do b=1,nV
        do c=1,nV

          do j=1,nO
            cWvovv(a,i,b,c) = cWvovv(a,i,b,c) + t1(j,a)*OOVV(j,i,b,c)
          end do

        end do
      end do
    end do
  end do

! OOOV block

  cWooov(:,:,:,:) = OOOV(:,:,:,:)

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do a=1,nO

          do b=1,nV
            cWooov(i,j,k,a) = cWooov(i,j,k,a) + t1(i,b)*OOVV(j,k,b,a)
          end do

        end do
      end do
    end do
  end do

! OOVV block

  cWoovv(:,:,:,:) = OOVV(:,:,:,:)

! VOOV block

  cWvoov(:,:,:,:) = VOOV(:,:,:,:)

  do a=1,nV
    do j=1,nO
      do i=1,nO
        do b=1,nV

          do e=1,nV
            cWvoov(a,j,i,b) = cWvoov(a,j,i,b) + t1(i,e)*VOOV(a,j,e,b)
          end do

          do m=1,nO
            cWvoov(a,j,i,b) = cWvoov(a,j,i,b) - t1(m,a)*OOOV(m,j,i,b)
          end do

          do m=1,nO
            do e=1,nO
              cWvoov(a,j,i,b) = cWvoov(a,j,i,b) - (t2(i,m,e,a) + t1(i,e)*t1(m,a))*OOVV(m,j,e,b)
            end do
          end do

        end do
      end do
    end do
  end do

! VVVO block

  cWvvvo(:,:,:,:) = VVVO(:,:,:,:)

  do a=1,nV
    do b=1,nV
      do c=1,nV
        do i=1,nO

          do j=1,nO
            do d=1,nV
              cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + t2(i,j,b,d)*VOVV(a,j,c,d)
              cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) - t2(i,j,a,d)*VOVV(b,j,c,d)
            end do
          end do

          do j=1,nO
            do k=1,nO
              cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + 0.5d0*tau(j,k,a,b)*VOOO(c,i,j,k)
            end do
          end do

          do j=1,nO
            cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + t2(i,j,a,b)*cFov(j,c)
          end do

          do d=1,nV
            cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + t1(i,d)*cWvvvv(a,b,c,d)
          end do

          do j=1,nO
            cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) - t1(j,a)*OVVO(j,b,c,i)
            cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + t1(j,b)*OVVO(j,a,c,i)

            do n=1,nO
              do d=1,nV
                cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) + t1(j,a)*t2(k,i,b,d)*OOVV(j,k,c,d)
                cWvvvo(a,b,c,i) = cWvvvo(a,b,c,i) - t1(j,b)*t2(k,i,a,d)*OOVV(j,k,c,d)
              end do
            end do

          end do

        end do
      end do
    end do
  end do

! OVOO block

  cWovoo(:,:,:,:) = OVOO(:,:,:,:)

  do i=1,nO
    do a=1,nV
      do j=1,nO
        do k=1,nO

          do l=1,nO
            do b=1,nV
              cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + t2(k,l,a,b)*OOOV(i,l,j,b)
              cWovoo(i,a,j,k) = cWovoo(i,a,j,k) - t2(k,l,a,b)*OOOV(j,l,i,b)
            end do
          end do

          do b=1,nV
            do c=1,nV
              cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + 0.5d0*tau(j,k,b,c)*OVVV(i,a,b,c)
            end do
          end do

          do b=1,nV
            cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + t2(j,k,a,b)*cFov(i,b)
          end do

          do l=1,nO
            cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + t1(l,a)*cWoooo(i,l,j,k)
          end do

          do b=1,nV
            cWovoo(i,a,j,k) = cWovoo(i,a,j,k) - t1(j,b)*OVVO(i,a,b,k)
            cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + t1(i,b)*OVVO(j,a,b,k) 

            do l=1,nO
              do c=1,nV
                cWovoo(i,a,j,k) = cWovoo(i,a,j,k) + t1(j,b)*t2(l,k,a,c)*OOVV(i,l,b,c)
                cWovoo(i,a,j,k) = cWovoo(i,a,j,k) - t1(i,b)*t2(l,k,a,c)*OOVV(j,l,b,c)
              end do
            end do

          end do

        end do
      end do
    end do
  end do

end subroutine 
