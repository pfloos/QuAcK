subroutine form_ring_r(nO,nV,OVVO,OOVV,t2,r2)

! Form residuals for ring CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OVVO(nO,nV,nV,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  r2(:,:,:,:) = 0d0

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV

          do k=1,nO
            do c=1,nV

              r2(i,j,a,b) = r2(i,j,a,b) + OVVO(i,c,a,k)*t2(k,j,c,b) + OVVO(j,c,b,k)*t2(i,k,a,c)

            end do
          end do

          do k=1,nO
            do l=1,nO
              do c=1,nV
                do d=1,nV

                  r2(i,j,a,b) = r2(i,j,a,b) + t2(i,k,a,c)*OOVV(k,l,c,d)*t2(l,j,d,b)

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_ring_r
