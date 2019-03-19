subroutine form_ub(nO,nV,OOVV,t1,ub)

! Form 1st term in (T) correction

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

  double precision,intent(in)   :: t1(nO,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: ub(nO,nO,nO,nV,nV,nV)

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do a=1,nV
          do b=1,nV
            do c=1,nV

              ub(i,j,k,a,b,c) = t1(i,a)*OOVV(j,k,b,c) &
                              + t1(i,b)*OOVV(j,k,c,a) &
                              + t1(i,c)*OOVV(j,k,a,b) &
                              + t1(j,a)*OOVV(k,i,b,c) & 
                              + t1(j,b)*OOVV(k,i,c,a) &
                              + t1(j,c)*OOVV(k,i,a,b) & 
                              + t1(k,a)*OOVV(i,j,b,c) &
                              + t1(k,b)*OOVV(i,j,c,a) & 
                              + t1(k,c)*OOVV(i,j,a,b)

            end do
          end do
        end do
      end do
    end do
  end do

end subroutine form_ub
