subroutine complex_vecout(n,v)

  implicit none

! Input variables
  integer,intent(in)            :: n
  complex*16,intent(in)         :: v(n)
! Local variables
  double precision,allocatable        :: v2(:,:)

  allocate(v2(n,2))
  write(*,*) 'First column real part, second imaginary part'
  v2(:,1) = real(v)
  v2(:,2) = aimag(v)
  call matout(n,2,v2)
  deallocate(v2)
end subroutine
