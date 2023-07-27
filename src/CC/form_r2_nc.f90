subroutine form_r2_nc(nO,nV,t1,t2,tau,delta_oovv,cFoo,cFov,cFvv,cWoooo,cWvvvv,cWovvo,OVOO,OOVV,OVVO,VVVO,r2)

! Form t2 residues in non-canonical CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: cFoo(nO,nO)
  double precision,intent(in)   :: cFov(nO,nV)
  double precision,intent(in)   :: cFvv(nV,nV)

  double precision,intent(in)   :: delta_oovv(nO,nO,nV,nV)

  double precision,intent(in)   :: cWoooo(nO,nO,nO,nO)
  double precision,intent(in)   :: cWvvvv(nV,nV,nV,nV)
  double precision,intent(in)   :: cWovvo(nO,nV,nV,nO)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

  double precision,intent(in)   :: OVOO(nO,nV,nO,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: VVVO(nV,nV,nV,nO)

! Local variables

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = OOVV(:,:,:,:)

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV

          do e=1,nV
            r2(i,j,a,b) = r2(i,j,a,b) + t2(i,j,a,e)*cFvv(b,e)
            r2(i,j,a,b) = r2(i,j,a,b) - t2(i,j,b,e)*cFvv(a,e)
          end do

          do m=1,nO
            do e=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,j,a,e)*t1(m,b)*cFov(m,e)
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(i,j,b,e)*t1(m,a)*cFov(m,e)
            end do
          end do

          do m=1,nO
            r2(i,j,a,b) = r2(i,j,a,b) - t2(i,m,a,b)*cFoo(m,j)
            r2(i,j,a,b) = r2(i,j,a,b) + t2(j,m,a,b)*cFoo(m,i)
          end do

          do m=1,nO
            do e=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,m,a,b)*t1(j,e)*cFov(m,e)
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(j,m,a,b)*t1(i,e)*cFov(m,e)
            end do
          end do

          do m=1,nO
            do n=1,nO
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*tau(m,n,a,b)*cWoooo(m,n,i,j)
            end do
          end do

          do e=1,nV
            do f=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*tau(i,j,e,f)*cWvvvv(a,b,e,f)
            end do
          end do

          do m=1,nO
            do e=1,nV
              r2(i,j,a,b) = r2(i,j,a,b)                                                 & 
                          + t2(i,m,a,e)*cWovvo(m,b,e,j) - t1(i,e)*t1(m,a)*OVVO(m,b,e,j) &
                          - t2(j,m,a,e)*cWovvo(m,b,e,i) + t1(j,e)*t1(m,a)*OVVO(m,b,e,i) &
                          - t2(i,m,b,e)*cWovvo(m,a,e,j) + t1(i,e)*t1(m,b)*OVVO(m,a,e,j) &
                          + t2(j,m,b,e)*cWovvo(m,a,e,i) - t1(j,e)*t1(m,b)*OVVO(m,a,e,i)
            end do
          end do

        do e=1,nV
          r2(i,j,a,b) = r2(i,j,a,b) + t1(i,e)*VVVO(a,b,e,j) - t1(j,e)*VVVO(a,b,e,i)
        end do

        do m=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - t1(m,a)*OVOO(m,b,i,j) + t1(m,b)*OVOO(m,a,i,j)
        end do

        end do
      end do
    end do
  end do

! Final expression of the t2 residue

  r2(:,:,:,:) = delta_oovv(:,:,:,:)*t2(:,:,:,:) - r2(:,:,:,:)

end subroutine 
