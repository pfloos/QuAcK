subroutine form_ubb(nO,nV,VVVO,VOOO,t2,ubb)

! Form 2nd term in (T) correction

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: VVVO(nV,nV,nV,nO)
  double precision,intent(in)   :: VOOO(nV,nO,nO,nO)

  double precision,intent(in)   :: t2(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,k,l,m
  integer                       :: a,b,c,d,e

! Output variables

  double precision,intent(out)  :: ubb(nO,nO,nO,nV,nV,nV)

  ubb(:,:,:,:,:,:) = 0d0

  do i=1,nO
    do j=1,nO
      do k=1,nO
        do a=1,nV
          do b=1,nV
            do c=1,nV

              do e=1,nV
                ubb(i,j,k,a,b,c) = ubb(i,j,k,a,b,c)          &
                                 + t2(i,j,a,e)*VVVO(b,c,e,k) & 
                                 + t2(i,j,b,e)*VVVO(c,a,e,k) & 
                                 + t2(i,j,c,e)*VVVO(a,b,e,k) & 
                                 + t2(k,i,a,e)*VVVO(b,c,e,j) & 
                                 + t2(k,i,b,e)*VVVO(c,a,e,j) & 
                                 + t2(k,i,c,e)*VVVO(a,b,e,j) & 
                                 + t2(j,k,a,e)*VVVO(b,c,e,i) & 
                                 + t2(j,k,b,e)*VVVO(c,a,e,i) & 
                                 + t2(j,k,c,e)*VVVO(a,b,e,i)
              end do
 
              do m=1,nO
                ubb(i,j,k,a,b,c) = ubb(i,j,k,a,b,c)          &
                                 + t2(i,m,a,b)*VOOO(c,m,j,k) & 
                                 + t2(i,m,b,c)*VOOO(a,m,j,k) & 
                                 + t2(i,m,c,a)*VOOO(b,m,j,k) & 
                                 + t2(j,m,a,b)*VOOO(c,m,k,i) & 
                                 + t2(j,m,b,c)*VOOO(a,m,k,i) & 
                                 + t2(j,m,c,a)*VOOO(b,m,k,i) & 
                                 + t2(k,m,a,b)*VOOO(c,m,i,j) & 
                                 + t2(k,m,b,c)*VOOO(a,m,i,j) & 
                                 + t2(k,m,c,a)*VOOO(b,m,i,j)
              end do
              
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine form_ubb
