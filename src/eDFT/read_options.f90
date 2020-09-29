subroutine read_options(nBas,method,x_rung,x_DFA,c_rung,c_DFA,SGn,nEns,wEns,aCC_w1,aCC_w2, & 
                        maxSCF,thresh,DIIS,max_diis,guess_type,ortho_type,doNcentered,occnum,Cx_choice)

! Read DFT options

  implicit none

  include 'parameters.h'

! Input variables
  integer,intent(in)           :: nBas

! Local variables

  integer                       :: iBas
  integer                       :: iEns
  integer                       :: iParam
  character(len=1)              :: answer
  double precision,allocatable  :: nEl(:)

! Output variables

  character(len=8),intent(out)  :: method
  integer,intent(out)           :: x_rung,c_rung
  character(len=12),intent(out) :: x_DFA, c_DFA
  integer,intent(out)           :: SGn
  integer,intent(out)           :: nEns
  logical,intent(out)           :: doNcentered
  double precision,intent(out)  :: wEns(maxEns)
  double precision,intent(out)  :: aCC_w1(3)
  double precision,intent(out)  :: aCC_w2(3)
  double precision,intent(out)  :: occnum(nBas,nspin,maxEns)

  integer,intent(out)           :: maxSCF
  double precision,intent(out)  :: thresh
  logical,intent(out)           :: DIIS
  integer,intent(out)           :: max_diis
  integer,intent(out)           :: guess_type
  integer,intent(out)           :: ortho_type
  integer,intent(out)           :: Cx_choice

! Open file with method specification

  open(unit=1,file='input/dft')

! Default values

  method  = 'GOK-RKS'
  x_rung  = 1
  c_rung  = 1
  x_DFA   = 'RS51'
  c_DFA   = 'RVWN5'
  SGn     = 0
  wEns(:) = 0d0

! Restricted or unrestricted calculation

  read(1,*)
  read(1,*) method

! EXCHANGE: read rung of Jacob's ladder

  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) x_rung,x_DFA

! CORRELATION: read rung of Jacob's ladder

  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
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
  
! Read occupation numbers for orbitals nO and nO+1

  occnum(:,:,:) = 0d0

  do iEns=1,nEns
    read(1,*)
    read(1,*) (occnum(iBas,1,iEns),iBas=1,nBas)
    read(1,*) (occnum(iBas,2,iEns),iBas=1,nBas)
  end do

  do iEns=1,nEns
    write(*,*) 
    write(*,*) '==============='
    write(*,*) 'State n.',iEns
    write(*,*) '==============='
    write(*,*)
    write(*,*) 'Spin-up   occupation numbers'
    write(*,*) (int(occnum(iBas,1,iEns)),iBas=1,nBas)
    write(*,*) 'Spin-down occupation numbers'
    write(*,*) (int(occnum(iBas,2,iEns)),iBas=1,nBas)
    write(*,*)
  end do
! Read ensemble weights for real physical (fractional number of electrons) ensemble (w1,w2)

  allocate(nEl(nEns))
  nEl(:) = 0d0
  do iEns=1,nEns
    do iBas=1,nBas
      nEl(iEns) = nEl(iEns) + occnum(iBas,1,iEns) + occnum(iBas,2,iEns)
    end do
  end do
  print*,'nEl'
  print*,nEl

  doNcentered = .false.

  read(1,*)
  read(1,*) (wEns(iEns),iEns=2,nEns)
  read(1,*)
  read(1,*) answer

  if(answer == 'T') doNcentered = .true.

  if (doNcentered) then

    wEns(1) = 1d0 - wEns(2) - wEns(3) 

  else

!   wEns(1) = 1d0 -  nEl(2)/nEl(1)*wEns(2) - nEl(3)/nEl(1)*wEns(3) 

    wEns(1) = 1d0 - wEns(2) - wEns(3)
    wEns(2) = nEl(1)/nEl(2)*wEns(2)
    wEns(3) = nEl(1)/nEl(3)*wEns(3)

  end if

  write(*,*)'----------------------------------------------------------'
  write(*,*)' Ensemble weights '
  write(*,*)'----------------------------------------------------------'
  call matout(nEns,1,wEns)
  write(*,*) 
  
! Read parameters for weight-dependent functional
  read(1,*)
  read(1,*) (aCC_w1(iParam),iParam=1,3)
  read(1,*) (aCC_w2(iParam),iParam=1,3)
! Read choice of exchange coefficient
  read(1,*)
  read(1,*) Cx_choice

  write(*,*)'----------------------------------------------------------'
  write(*,*)' parameters for w1-dependant exchange functional coefficient '
  write(*,*)'----------------------------------------------------------'
  call matout(3,1,aCC_w1)
  write(*,*)

  write(*,*)'----------------------------------------------------------'
  write(*,*)' parameters for w2-dependant exchange functional coefficient '
  write(*,*)'----------------------------------------------------------'
  call matout(3,1,aCC_w2)
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
