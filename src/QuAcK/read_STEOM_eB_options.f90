subroutine read_STEOM_eB_options(working_dir,TDA,singlet,triplet,Davidson,nRoots,tolerance,model)

  use, intrinsic :: iso_fortran_env, only: iostat_end
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  character(*),intent(in) :: working_dir
  logical,intent(out) :: TDA,singlet,triplet,Davidson
  character(*),intent(out) :: model
  integer,intent(out) :: nRoots
  double precision,intent(out) :: tolerance
  character(512) :: line
  character(1) :: flags(4)
  integer :: unit,status

  TDA=.false.
  singlet=.true.
  triplet=.true.
  Davidson=.false.
  nRoots=5
  tolerance=1d-8
  model='COT2'

  open(newunit=unit,file=trim(working_dir)//'/input/options',status='old',action='read',iostat=status)
  if(status /= 0) error stop 'Cannot open options for STEOM-eB'
  do
    read(unit,'(A)',iostat=status) line
    if(status == iostat_end) exit
    if(status /= 0) error stop 'Error reading STEOM-eB options'
    if(index(adjustl(line),'# STEOM-eB:') /= 1) cycle

    read(unit,'(A)',iostat=status) line
    if(status /= 0) error stop 'Missing STEOM-eB option values'
    flags=' '
    read(line,*,iostat=status) flags(1:3)
    if(status /= 0 .or. any(flags(1:3) /= 'T' .and. flags(1:3) /= 'F')) &
      error stop 'Expected T or F for STEOM-eB TDA, singlet and triplet'

    TDA=flags(1) == 'T'
    singlet=flags(2) == 'T'
    triplet=flags(3) == 'T'

    read(line,*,iostat=status) flags
    if(status == iostat_end) exit
    if(status /= 0 .or. (flags(4) /= 'T' .and. flags(4) /= 'F')) &
      error stop 'Expected T or F for STEOM-eB Davidson'
    Davidson=flags(4) == 'T'
    read(line,*,iostat=status) flags,nRoots,tolerance
    if(status /= 0) error stop 'Expected STEOM-eB Davidson root count and threshold'
    read(line,*,iostat=status) flags,nRoots,tolerance,model
    if(status /= 0 .and. status /= iostat_end) error stop 'Error reading the STEOM-eB model'
    exit
  end do
  close(unit)

  if(model /= 'COT2' .and. model /= 'legacy' .and. model /= 'osCOT') &
    error stop 'Expected COT2, legacy or osCOT for the STEOM-eB model'
  if(nRoots <= 0 .or. .not.ieee_is_finite(tolerance) .or. tolerance <= 0d0) &
    error stop 'STEOM-eB Davidson requires a positive root count and finite positive threshold'
  if(.not.(singlet.or.triplet)) error stop 'STEOM-eB requires at least one spin manifold'

end subroutine read_STEOM_eB_options
