subroutine form_ladder_r(nC,nO,nV,nR,OOOO,OOVV,VVVV,t2,r2)

! Form residuals for ladder CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: y(:,:,:,:)

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = 0d0

  allocate(y(nO,nO,nO,nO))

  y(:,:,:,:) = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do k=nC+1,nO
        do l=nC+1,nO
          do c=1,nV-nR
            do d=1,nV-nR
              y(i,j,k,l) = y(i,j,k,l) + t2(i,j,c,d)*OOVV(k,l,c,d)
            end do
          end do
        end do
      end do
    end do
  end do

  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR

          do k=nC+1,nO
            do l=nC+1,nO
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(k,l,a,b)*OOOO(i,j,k,l)
            end do
          end do

          do c=1,nV-nR
            do d=1,nV-nR
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*VVVV(c,d,a,b)*t2(i,j,c,d) 
            end do
          end do

          do k=nC+1,nO
            do l=nC+1,nO
              r2(i,j,a,b) = r2(i,j,a,b) + 0.25d0*y(i,j,k,l)*t2(k,l,a,b)
            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_ladder_r
