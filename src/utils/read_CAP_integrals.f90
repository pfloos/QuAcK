subroutine read_CAP_integrals(nBas,W)

! Read one- and two-electron integrals from files

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: debug
  integer                       :: mu,nu
  double precision              :: wxyz

! Output variables

  double precision,intent(out)  :: W(nBas,nBas)

! Open file with integrals

  debug = .false.
  open(unit=31,file='int/CAP.dat')

! Read CAP integrals

  W(:,:) = 0d0
  do 
    read(31,*,end=31) mu,nu,wxyz
    W(mu,nu) = wxyz
    W(nu,mu) = wxyz
  end do
  31 close(unit=31)

! Print results
  if(debug) then
    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'CAP integrals'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,W)
    write(*,*)
  end if

end subroutine 
