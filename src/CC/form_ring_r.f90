subroutine form_ring_r(nC,nO,nV,nR,OVVO,VOOV,OOVV,t2,r2)

! Form residuals for ring CCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: VOOV(nV,nO,nO,nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: y(:,:,:,:)

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = 0d0

  allocate(y(nO,nV,nO,nV))

  y(:,:,:,:) = 0d0

  do i=nC+1,nO
    do a=1,nV-nR
      do l=nC+1,nO
        do d=1,nV-nR
          do k=nC+1,nO
            do c=1,nV-nR
              y(i,a,l,d) = y(i,a,l,d) + t2(i,k,a,c)*OOVV(k,l,c,d)
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

              r2(i,j,a,b) = r2(i,j,a,b) + VOOV(a,k,i,c)*t2(k,j,c,b) + OVVO(k,b,c,j)*t2(i,k,a,c)

            end do
          end do

          do l=nC+1,nO
            do d=1,nV-nR

              r2(i,j,a,b) = r2(i,j,a,b) + y(i,a,l,d)*t2(l,j,d,b)

            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_ring_r
