subroutine read_faddeev_options(working_dir,reference,pairs,coupling,static3,external_gw,davidson,nroots)
  implicit none
  character(*),intent(in) :: working_dir
  integer,intent(out) :: reference,pairs,coupling,nroots
  logical,intent(out) :: static3,external_gw,davidson
  character(512) :: line
  character(1) :: flags(3)
  integer :: unit,status

  reference=0; pairs=0; coupling=3; nroots=1
  static3=.true.; external_gw=.false.; davidson=.true.
  open(newunit=unit,file=trim(working_dir)//'/input/options',status='old',action='read',iostat=status)
  if(status/=0) error stop 'Cannot open Faddeev options'
  do
    read(unit,'(A)',iostat=status) line
    if(status<0) exit
    if(status/=0) error stop 'Error reading Faddeev options'
    if(index(adjustl(line),'# Faddeev:')/=1) cycle
    read(unit,'(A)',iostat=status) line
    if(status/=0) error stop 'Missing Faddeev option values'
    flags=' '
    read(line,*,iostat=status) reference,pairs,coupling,flags,nroots
    if(status/=0.or.any(flags/='T'.and.flags/='F')) error stop 'Invalid Faddeev option row'
    static3=flags(1)=='T'; external_gw=flags(2)=='T'; davidson=flags(3)=='T'
    exit
  end do
  close(unit)
  if(reference<0.or.reference>1.or.pairs<0.or.pairs>2.or.coupling<2.or.coupling>5.or.nroots<1) &
    error stop 'Faddeev option out of range'
  if(pairs==1.and.reference/=1) error stop 'BSE pair interactions require a GW reference'
end subroutine
