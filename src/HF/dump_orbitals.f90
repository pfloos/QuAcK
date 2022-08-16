subroutine dump_orbitals(nBas,c)

! Write orbitals in a file

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  integer                       :: mu,p

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

! Open file with orbitals

  open(unit=12,file='int/MO.dat')

! Dump orbitals

  do mu=1,nBas
    do p=1,nBas
      if(abs(c(mu,p)) > 1d-15) write(12,'(I6,I6,F20.15)') mu,p,c(mu,p)
    enddo
  enddo

! Close file

  12 close(unit=12)

end subroutine dump_orbitals
