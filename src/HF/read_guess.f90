subroutine read_guess(nBas,c)

! Read the orbitals from a file

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: debug
  integer                       :: mu,p
  double precision              :: coeff

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

! Open file with orbitals

  debug = .false.

  open(unit=12,file='int/MO.dat')

! Read coefficients

  c(:,:) = 0d0
  do 
    read(12,*,end=12) mu,p,coeff
    c(mu,p) = coeff
    c(mu,p) = coeff
  enddo

! Close file

  12 close(unit=12)

! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Orbitals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,c)
    write(*,*)
  endif

end subroutine read_guess
