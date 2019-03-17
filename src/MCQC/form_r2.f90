subroutine form_r2(nO,nV,OOVV,OVOO,OVVV,OVVO,gvv,goo,aoooo,bvvvv,hovvo,t1,t2,tau,r2)

! Form tau in CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVOO(nO,nV,nO,nO)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)

  double precision,intent(in)   :: gvv(nV,nV)
  double precision,intent(in)   :: goo(nO,nO)
  double precision,intent(in)   :: aoooo(nO,nO,nO,nO)
  double precision,intent(in)   :: bvvvv(nV,nV,nV,nV)
  double precision,intent(in)   :: hovvo(nO,nV,nV,nO)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = OOVV(:,:,:,:)

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV

          do k=1,nO
            do l=1,nO
              r2(i,j,a,b) = r2(i,j,a,b) + aoooo(i,j,k,l)*tau(k,l,a,b)
            end do
          end do

          do c=1,nV
            do d=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + bvvvv(c,d,a,b)*tau(i,j,c,d)
            end do
          end do

          do c=1,nV
            r2(i,j,a,b) = r2(i,j,a,b) + gvv(c,a)*t2(i,j,c,b)
          end do

          do k=1,nO
            r2(i,j,a,b) = r2(i,j,a,b) + OVOO(k,a,i,j)*t1(k,b)
          end do

          do c=1,nV
            r2(i,j,a,b) = r2(i,j,a,b) - gvv(c,b)*t2(i,j,c,a)
          end do

          do k=1,nO
            r2(i,j,a,b) = r2(i,j,a,b) - OVOO(k,b,i,j)*t1(k,a)
          end do

          do k=1,nO
            r2(i,j,a,b) = r2(i,j,a,b) - goo(i,k)*t2(k,j,a,b)
          end do

          do c=1,nV
            r2(i,j,a,b) = r2(i,j,a,b) + OVVV(j,c,b,a)*t1(i,c)
          end do

          do k=1,nO
            r2(i,j,a,b) = r2(i,j,a,b) + goo(j,k)*t2(k,i,a,b)
          end do

          do c=1,nV
            r2(i,j,a,b) = r2(i,j,a,b) - OVVV(i,c,b,a)*t1(j,c)
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + hovvo(i,c,a,k)*t2(j,k,b,c)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - OVVO(i,c,a,k)*t1(j,c)*t1(k,b)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - hovvo(j,c,a,k)*t2(i,k,b,c)
            end do
          end do

          do k=1,nO
            do c=1,nV
               r2(i,j,a,b) = r2(i,j,a,b) + OVVO(j,c,a,k)*t1(i,c)*t1(k,b)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - hovvo(i,c,b,k)*t2(j,k,a,c)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + OVVO(i,c,b,k)*t1(j,c)*t1(k,a)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + hovvo(j,c,b,k)*t2(i,k,a,c)
            end do
          end do

          do k=1,nO
            do c=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) - OVVO(j,c,b,k)*t1(i,c)*t1(k,a)
            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_r2
