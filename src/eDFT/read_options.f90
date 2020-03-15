subroutine read_options(method,x_rung,x_DFA,c_rung,c_DFA,SGn,nEns,wEns,maxSCF,thresh,DIIS,max_diis,guess_type,ortho_type)

! Read DFT options

  implicit none

  include 'parameters.h'

! Local variables

  integer                       :: I

! Output variables

  character(len=7),intent(out)  :: method
  integer,intent(out)           :: x_rung,c_rung
  character(len=12),intent(out)  :: x_DFA, c_DFA
  integer,intent(out)           :: SGn
  integer,intent(out)           :: nEns
  double precision,intent(out)  :: wEns(maxEns)

  integer,intent(out)           :: maxSCF
  double precision,intent(out)  :: thresh
  logical,intent(out)           :: DIIS
  integer,intent(out)           :: max_diis
  integer,intent(out)           :: guess_type
  integer,intent(out)           :: ortho_type

! Local variables

  character(len=1)              :: answer

! Open file with method specification

  open(unit=1,file='input/dft')

! Default values

  method  = 'GOK-RKS'
  x_rung  = 1
  c_rung  = 1
  x_DFA   = 'S51'
  c_DFA   = 'W38'
  SGn     = 0
  wEns(:) = 0d0

! Restricted or unrestricted calculation

  read(1,*)
  read(1,*) method

! EXCHANGE: read rung of Jacob's ladder

  read(1,*)
  read(1,*) x_rung,x_DFA

! CORRELATION: read rung of Jacob's ladder

  read(1,*)
  read(1,*) c_rung,c_DFA

! Read SG-n grid

  read(1,*)
  read(1,*) SGn

! Read number of states in ensemble

  read(1,*)
  read(1,*) nEns

  if(nEns.gt.maxEns) then
    write(*,*) ' Number of states in ensemble too big!! ' 
    stop
  endif

  write(*,*)'----------------------------------------------------------'
  write(*,'(A33,I3)')'  Number of states in ensemble = ',nEns
  write(*,*)'----------------------------------------------------------'
  write(*,*) 
  
! Read ensemble weights
  read(1,*)
  read(1,*) (wEns(I),I=2,nEns)
  wEns(1) = 1d0 - sum(wEns)

  write(*,*)'----------------------------------------------------------'
  write(*,*)' Ensemble weights '
  write(*,*)'----------------------------------------------------------'
  call matout(nEns,1,wEns)
  write(*,*) 

! Read KS options

  maxSCF     = 64
  thresh     = 1d-6
  DIIS       = .false.
  max_diis   = 5
  guess_type = 1
  ortho_type = 1

  read(1,*)
  read(1,*) maxSCF,thresh,answer,max_diis,guess_type,ortho_type

  if(answer == 'T') DIIS = .true.

  if(.not.DIIS) max_diis = 1

! Close file with options

  close(unit=1)

end subroutine read_options
