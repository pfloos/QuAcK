subroutine GWiw_read_wweight(nfreqs,wweight)

! Read the imaginary-frequency quadrature weights from ./grids/wweight.txt (greenX minimax grid,
! same file and format as read_scGX_grids).

  implicit none

! Input variables

  integer,intent(in)            :: nfreqs

! Output variables

  double precision,intent(out)  :: wweight(nfreqs)

! Local variables

  logical                       :: file_exists
  integer                       :: iunit,ivar,ifreq

  inquire(file='./grids/wweight.txt',exist=file_exists)
  if(.not.file_exists) then
    write(*,*) ' Error! Could not find ./grids/wweight.txt (needed for EcGM from imaginary frequencies)'
    error stop
  end if

  open(newunit=iunit,form='formatted',file='./grids/wweight.txt',status='old')
  read(iunit,*) ivar
  if(ivar /= nfreqs) then
    write(*,*) ' Error! nfreqs in the options file does not match ./grids/wweight.txt: ',nfreqs,ivar
    error stop
  end if
  do ifreq=1,nfreqs
    read(iunit,*) wweight(ifreq)
  end do
  close(iunit)

end subroutine
