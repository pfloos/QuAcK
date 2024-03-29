subroutine form_h(nC,nO,nV,nR,eO,eV,OOVV,t1,tau,hvv,hoo,hvo)

! Scuseria Eqs. (5), (6) and (7)

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR

  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: hvv(nV,nV)
  double precision,intent(out)  :: hoo(nO,nO)
  double precision,intent(out)  :: hvo(nV,nO)

  hvv(:,:) = 0d0

  do b=1,nV-nR
    hvv(b,b) = eV(b)
    do a=1,nV-nR
      do j=nC+1,nO
        do k=nC+1,nO
          do c=1,nV-nR

            hvv(b,a) = hvv(b,a) - OOVV(j,k,b,c)*tau(j,k,a,c)

          end do
        end do
      end do
    end do
  end do

  hoo(:,:) = 0d0

  do i=nC+1,nO
    hoo(i,i) = eO(i)
    do j=nC+1,nO
      do k=nC+1,nO
        do b=1,nV-nR
          do c=1,nV-nR

            hoo(i,j) = hoo(i,j) + OOVV(j,k,b,c)*tau(i,k,b,c)

          end do
        end do
      end do
    end do
  end do

  hvo(:,:) = 0d0

  do b=1,nV-nR
    do j=nC+1,nO
      do k=nC+1,nO
        do c=1,nV-nR

          hvo(b,j) = hvo(b,j) + OOVV(j,k,b,c)*t1(k,c)

        end do
      end do
    end do
  end do

! print*,'hvv',hvv

end subroutine 
