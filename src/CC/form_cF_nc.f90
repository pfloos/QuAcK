subroutine form_cF_nc(nO,nV,t1,taus,Foo,Fov,Fvv,OOOV,OOVV,OVVV,cFoo,cFov,cFvv)

! Compute F terms in CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: taus(nO,nO,nV,nV)

  double precision,intent(in)   :: Foo(nO,nO)
  double precision,intent(in)   :: Fov(nO,nV)
  double precision,intent(in)   :: Fvv(nV,nV)

  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)


! Local variables

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision,external     :: Kronecker_Delta

! Output variables

  double precision,intent(out)  :: cFoo(nO,nO)
  double precision,intent(out)  :: cFov(nO,nV)
  double precision,intent(out)  :: cFvv(nV,nV)

! Occupied-occupied block

  do m=1,nO
    do i=1,nO

      cFoo(m,i) = (1d0 - Kronecker_delta(m,i))*Foo(m,i) 

      do e=1,nV
        cFoo(m,i) = cFoo(m,i) + 0.5d0*t1(i,e)*Fov(m,e)
      end do

      do e=1,nV
        do n=1,nO
          cFoo(m,i) = cFoo(m,i) + t1(n,e)*OOOV(m,n,i,e)
        end do
      end do

      do e=1,nV
        do n=1,nO
          do f=1,nV
            cFoo(m,i) = cFoo(m,i) + 0.5d0*taus(i,n,e,f)*OOVV(m,n,e,f)
          end do
        end do
      end do

    end do
  end do

! Occupied-virtual block

  cFov(:,:) = Fov(:,:)

  do m=1,nO
    do e=1,nV

      do n=1,nO
        do f=1,nV
          cFov(m,e) = cFov(m,e) + t1(n,f)*OOVV(m,n,e,f)
        end do
      end do

    end do
  end do

! Virtual-virtual block

  do a=1,nV
    do e=1,nV

      cFvv(a,e) = (1d0 - Kronecker_delta(a,e))*Fvv(a,e) 
 
      do m=1,nO
        cFvv(a,e) = cFvv(a,e) - 0.5d0*t1(m,a)*Fov(m,e)
      end do

      do m=1,nO
        do f=1,nV
        cFvv(a,e) = cFvv(a,e) + t1(m,f)*OVVV(m,a,f,e)
        end do
      end do

      do m=1,nO
        do n=1,nO
          do f=1,nV
           cFvv(a,e) = cFvv(a,e) - 0.5d0*taus(m,n,a,f)*OOVV(m,n,e,f)
          end do
        end do
      end do

    end do
  end do

end subroutine 
