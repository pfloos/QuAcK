subroutine form_g(nC,nO,nV,nR,hvv,hoo,VOVV,OOOV,t1,gvv,goo)

! Scuseria Eqs. (9), (10)

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR

  double precision,intent(in)   :: hvv(nV,nV)
  double precision,intent(in)   :: hoo(nO,nO)

  double precision,intent(in)   :: VOVV(nV,nO,nV,nV)
  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)

  double precision,intent(in)   :: t1(nO,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: gvv(nV,nV)
  double precision,intent(out)  :: goo(nO,nO)

  gvv(:,:) = hvv(:,:)

  do c=1,nV-nR
    do a=1,nV-nR
      do k=nC+1,nO
        do d=1,nV-nR
          gvv(c,a) = gvv(c,a) + VOVV(a,k,c,d)*t1(k,d)
        end do
      end do
    end do
  end do

  goo(:,:) = hoo(:,:)

  do i=nC+1,nO
    do k=nC+1,nO
      do l=nC+1,nO
        do c=1,nV-nR
          goo(i,k) = goo(i,k) + OOOV(k,l,i,c)*t1(l,c)
        end do
      end do
    end do
  end do

end subroutine 
