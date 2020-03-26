subroutine form_r1(nC,nO,nV,nR,OVVO,OVVV,OOOV,hvv,hoo,hvo,t1,t2,tau,r1)

! Form tau in CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR

  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)
  double precision,intent(in)   :: OOOV(nO,nO,nO,nV)


  double precision,intent(in)   :: hvv(nV,nV)
  double precision,intent(in)   :: hoo(nO,nO)
  double precision,intent(in)   :: hvo(nV,nO)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: r1(nO,nV)

  r1(:,:) = 0d0

  do i=nC+1,nO
      do a=1,nV-nR

        do b=1,nV-nR
          r1(i,a) = r1(i,a) + hvv(b,a)*t1(i,b)
        end do

        do j=nC+1,nO
          r1(i,a) = r1(i,a) - hoo(i,j)*t1(j,a)
        end do

        do j=nC+1,nO
          do b=1,nV-nR
            r1(i,a) = r1(i,a) + hvo(b,j)*(t2(i,j,a,b) + t1(i,b)*t1(j,a))
          end do
        end do

        do j=nC+1,nO
          do b=1,nV-nR
            r1(i,a) = r1(i,a) + OVVO(i,b,a,j)*t1(j,b)
          end do
        end do

        do j=nC+1,nO
          do b=1,nV-nR
            do c=1,nV-nR
              r1(i,a) = r1(i,a) - OVVV(j,a,b,c)*tau(i,j,b,c)
            end do
          end do
        end do

        do j=nC+1,nO
          do k=nC+1,nO
            do b=1,nV-nR
              r1(i,a) = r1(i,a) - OOOV(j,k,i,b)*tau(j,k,a,b)
            end do
          end do
        end do

    end do
  end do

end subroutine form_r1
