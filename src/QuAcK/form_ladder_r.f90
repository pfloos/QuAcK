subroutine form_ladder_r(nO,nV,OOOO,OOVV,VVVV,t2,r2)

! Form residuals for ladder CCD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)

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
            do l=1,nO
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(k,l,a,b)*OOOO(i,j,k,l)
            end do
          end do

          do c=1,nV
            do d=1,nV
              r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*VVVV(c,d,a,b)*t2(i,j,c,d) 
            end do
          end do

          do k=1,nO
            do l=1,nO
              do c=1,nV
                do d=1,nV
                  r2(i,j,a,b) = r2(i,j,a,b) + 0.25d0*t2(i,j,c,d)*OOVV(k,l,c,d)*t2(k,l,a,b)

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_ladder_r
