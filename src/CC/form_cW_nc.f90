subroutine form_cW_nc(nO,nV,t1,t2,tau,OOOO,OOOV,OOVO,OOVV,OVVO,OVVV,VOVV,VVVV,cWoooo,cWovvo,cWvvvv)

! Compute W terms in CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)
  double precision,intent(in)   :: OOVO(nO,nO,nV,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)
  double precision,intent(in)   :: VOVV(nV,nO,nV,nV)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)

! Local variables

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision,external     :: Kronecker_Delta

! Output variables

  double precision,intent(out)  :: cWoooo(nO,nO,nO,nO)
  double precision,intent(out)  :: cWovvo(nO,nV,nV,nO)
  double precision,intent(out)  :: cWvvvv(nV,nV,nV,nV)

! OOOO block  

  cWoooo(:,:,:,:) = OOOO(:,:,:,:)

  do m=1,nO
    do n=1,nO
      do i=1,nO
        do j=1,nO

          do e=1,nV
            cWoooo(m,n,i,j) = cWoooo(m,n,i,j) + t1(j,e)*OOOV(m,n,i,e) - t1(i,e)*OOOV(m,n,j,e)
          end do

          do e=1,nV
            do f=1,nV
            cWoooo(m,n,i,j) = cWoooo(m,n,i,j) + 0.25d0*tau(i,j,e,f)*OOVV(m,n,e,f)
            end do
          end do

        end do
      end do
    end do
  end do

! OVVO block

  cWovvo(:,:,:,:) = OVVO(:,:,:,:)

  do m=1,nO
    do b=1,nV
      do e=1,nV
        do j=1,nO
         
          do f=1,nV
            cWovvo(m,b,e,j) = cWovvo(m,b,e,j) + t1(j,f)*OVVV(m,b,e,f)
          end do

          do n=1,nO
            cWovvo(m,b,e,j) = cWovvo(m,b,e,j) - t1(n,b)*OOVO(m,n,e,j)
          end do

          do n=1,nO
            do f=1,nV
              cWovvo(m,b,e,j) = cWovvo(m,b,e,j) &
                              - ( 0.5d0*t2(j,n,f,b) + t1(j,f)*t1(n,b) )*OOVV(m,n,e,f)
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
      do e=1,nV
        do f=1,nV

          do m=1,nO
            cWvvvv(a,b,e,f) = cWvvvv(a,b,e,f) - t1(m,b)*VOVV(a,m,e,f) + t1(m,a)*VOVV(b,m,e,f)
          end do

          do m=1,nO
            do n=1,nO
              cWvvvv(a,b,e,f) = cWvvvv(a,b,e,f) + 0.25d0*tau(m,n,a,b)*OOVV(m,n,e,f)
            end do
          end do

        end do
      end do
    end do
  end do


end subroutine form_cW_nc
