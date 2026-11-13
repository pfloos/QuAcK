subroutine GWiw_read_wcoord(nfreqs,wcoord)

! Read the imaginary-frequency sampling points from ./grids/wcoord.txt (greenX minimax grid, same file
! and format as read_scGX_grids). Only the frequencies are needed by the contour-deformation SRG scheme.

  implicit none

! Input variables

  integer,intent(in)            :: nfreqs

! Output variables

  double precision,intent(out)  :: wcoord(nfreqs)

! Local variables

  logical                       :: file_exists
  integer                       :: iunit,ivar,ifreq

  inquire(file='./grids/wcoord.txt',exist=file_exists)
  if(.not.file_exists) then
    write(*,*) ' Error! Could not find ./grids/wcoord.txt (needed for the imaginary-axis SRG self-energy)'
    error stop
  end if

  open(newunit=iunit,form='formatted',file='./grids/wcoord.txt',status='old')
  read(iunit,*) ivar
  if(ivar /= nfreqs) then
    write(*,*) ' Error! nfreqs in the options file does not match ./grids/wcoord.txt: ',nfreqs,ivar
    error stop
  end if
  do ifreq=1,nfreqs
    read(iunit,*) wcoord(ifreq)
  end do
  close(iunit)

end subroutine
