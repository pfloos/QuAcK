subroutine form_crossed_ring_r(nC,nO,nV,nR,OVOV,OOVV,t2,r2)

! Form residuals for crossed-ring CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: y(:,:,:,:)

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = 0d0

  allocate(y(nV,nO,nO,nV))

  y(:,:,:,:) = 0d0

  do i=nC+1,nO
    do b=1,nV-nR
      do l=nC+1,nO
        do d=1,nV-nR
          do k=nC+1,nO
            do c=1,nV-nR
              y(b,i,l,d) = y(b,i,l,d) + t2(i,k,c,b)*OOVV(k,l,c,d)
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
            do c=1,nV-nR

              r2(i,j,a,b) = r2(i,j,a,b) - OVOV(k,b,i,c)*t2(k,j,a,c) - OVOV(k,a,j,c)*t2(i,k,c,b)

            end do
          end do

          do l=nC+1,nO
            do d=1,nV-nR

              r2(i,j,a,b) = r2(i,j,a,b) - y(b,i,l,d)*t2(l,j,a,d)

            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_crossed_ring_r
