subroutine read_COT_options(working_dir,TDA_W,useDavidson,TDA_COT,davidsonThreshold)

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  character(*),intent(in) :: working_dir
  logical,intent(in) :: TDA_W
  logical,intent(out) :: useDavidson,TDA_COT
  double precision,intent(out) :: davidsonThreshold
  character(512) :: line
  character(1) :: flags(2)
  integer :: unit,status

  useDavidson=.true.
  TDA_COT=TDA_W
  davidsonThreshold=1d-10
  open(newunit=unit,file=trim(working_dir)//'/input/options',status='old',action='read',iostat=status)
  if(status /= 0) error stop 'Cannot open COT options'
  do
    read(unit,'(A)',iostat=status) line
    if(status < 0) exit
    if(status /= 0) error stop 'Error reading COT options'
    if(index(adjustl(line),'# COT:') /= 1) cycle
    read(unit,'(A)',iostat=status) line
    if(status /= 0) error stop 'Missing COT option values'
    flags=' '
    read(line,*,iostat=status) flags
    if(status /= 0 .or. any(flags /= 'T' .and. flags /= 'F')) error stop 'Expected T or F for both COT options'
    read(line,*,iostat=status) flags,davidsonThreshold
    if(status /= 0 .and. .not.is_iostat_end(status)) error stop 'Invalid COT Davidson threshold'
    if(.not.ieee_is_finite(davidsonThreshold) .or. davidsonThreshold <= 0d0) &
      error stop 'COT Davidson threshold must be finite and positive'
    useDavidson=flags(1) == 'T'
    TDA_COT=flags(2) == 'T'
    exit
  end do
  close(unit)

end subroutine read_COT_options
