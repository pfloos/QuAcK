subroutine form_abh(nO,nV,OOOO,OVOO,OOVV,VVVV,VOVV,OVVO,OVVV,t1,tau,aoooo,bvvvv,hovvo)

! Scuseria Eqs. (11),(12) and (13)

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: OVOO(nO,nV,nO,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)
  double precision,intent(in)   :: VOVV(nV,nO,nV,nV)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: aoooo(nO,nO,nO,nO)
  double precision,intent(out)  :: bvvvv(nV,nV,nV,nV)
  double precision,intent(out)  :: hovvo(nO,nV,nV,nO)

  aoooo(:,:,:,:) = OOOO(:,:,:,:)

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do l=1,nO

          do c=1,nV
            aoooo(i,j,k,l) = aoooo(i,j,k,l) + OVOO(i,c,k,l)*t1(j,c)
          end do

          do c=1,nV
            aoooo(i,j,k,l) = aoooo(i,j,k,l) - OVOO(j,c,k,l)*t1(i,c)
          end do

          do c=1,nV
            do d=1,nV
              aoooo(i,j,k,l) = aoooo(i,j,k,l) + OOVV(k,l,c,d)*tau(i,j,c,d)
            end do
          end do

        end do
      end do
    end do
  end do

  bvvvv(:,:,:,:) = VVVV(:,:,:,:)

  do c=1,nV 
    do d=1,nV
      do a=1,nV
        do b=1,nV
          
          do k=1,nO
            bvvvv(c,d,a,b) = bvvvv(c,d,a,b) - VOVV(a,k,c,d)*t1(k,b)
          end do

          do k=1,nO
            bvvvv(c,d,a,b) = bvvvv(c,d,a,b) + VOVV(b,k,c,d)*t1(k,a)
          end do

        end do
      end do
    end do
  end do

  hovvo(:,:,:,:) = OVVO(:,:,:,:)

  do i=1,nO 
    do c=1,nV
      do a=1,nV
        do k=1,nO

          do l=1,nO
            hovvo(i,c,a,k) = hovvo(i,c,a,k) - OVOO(i,c,l,k)*t1(l,a)
          end do

          do d=1,nV
            hovvo(i,c,a,k) = hovvo(i,c,a,k) + OVVV(k,a,c,d)*t1(i,d)
          end do

          do l=1,nO
            do d=1,nV
              hovvo(i,c,a,k) = hovvo(i,c,a,k) - OOVV(k,l,c,d)*tau(i,l,d,a)
            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_abh
