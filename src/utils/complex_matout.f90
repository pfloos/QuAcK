subroutine complex_matout(n,m,A)

  implicit none

! Input variables

  integer,intent(in)            :: n,m
  complex*16,intent(in)         :: A(n,m)

  write( *,*) 'Real part:'
  call matout(n,m,real(A))
  write (*,*) 'Imaginary part:'
  call matout(n,m,aimag(A))

end subroutine
