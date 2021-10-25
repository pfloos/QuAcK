subroutine read_options_dft(nBas,method,x_rung,x_DFA,c_rung,c_DFA,SGn,nEns,wEns,aCC_w1,aCC_w2, & 
                            doNcentered,occnum,Cx_choice)

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
  character(len=12)             :: x_func,c_func

! Output variables

  character(len=8),intent(out)  :: method
  integer,intent(out)           :: x_rung,c_rung
  integer,intent(out)           :: x_DFA,c_DFA
  integer,intent(out)           :: SGn
  integer,intent(out)           :: nEns
  logical,intent(out)           :: doNcentered
  double precision,intent(out)  :: wEns(maxEns)
  double precision,intent(out)  :: aCC_w1(3)
  double precision,intent(out)  :: aCC_w2(3)
  double precision,intent(out)  :: occnum(nBas,nspin,maxEns)

  integer,intent(out)           :: Cx_choice

! Open file with method specification

  open(unit=1,file='input/dft')

! Default values

  method  = 'eDFT-UKS'
  x_rung  = 1
  c_rung  = 1
  x_DFA   = 1
  c_DFA   = 1
  SGn     = 0
  wEns(:) = 0d0

! Restricted or unrestricted calculation

  read(1,*)
  read(1,*) method

!---------------------------------------!
! EXCHANGE: read rung of Jacob's ladder !
!---------------------------------------!

  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) x_rung,x_func

  select case (x_rung) ! exchange functionals
 
    case (0) ! Hartree
   
    select case (x_func)

      case ('H')
     
        x_DFA = 1
     
      case default
     
        call print_warning('!!! Hartree exchange functional not available !!!')
        stop

    end select

    case (1) ! LDA
   
    select case (x_func)

      case ('S51')
     
        x_DFA = 1
     
      case ('CC-S51')
     
        x_DFA = 2
     
      case default
     
        call print_warning('!!! LDA exchange functional not available !!!')
        stop

    end select

    case (2) ! GGA
   
      select case (x_func)
 
        case ('G96')
  
          x_DFA = 1
  
        case ('B88')
  
          x_DFA = 2
  
        case ('PBE')
  
          x_DFA = 3
  
        case default
 
          call print_warning('!!! GGA exchange functional not available !!!')
          stop

      end select

    case (3) ! MGGA
   
    select case (x_func)

      case default

        call print_warning('!!! MGGA exchange functional not available !!!')
        stop

    end select

    case (4) ! Hybrid
   
      select case (x_func)
 
        case ('HF')
  
          x_DFA = 1
  
        case ('B3LYP')
  
          x_DFA = 2
  
        case ('BHHLYP')
 
          x_DFA = 3
  
        case ('PBE')
  
          x_DFA = 4
  
        case default
 
          call print_warning('!!! Hybrid exchange functional not available !!!')
          stop
 
      end select

    case default

      call print_warning('!!! Exchange rung not available !!!')
      stop

  end select

! Select rung for exchange 

  write(*,*)
  write(*,*) '*******************************************************************'
  write(*,*) '*                        Exchange rung                            *'
  write(*,*) '*******************************************************************'

  call select_rung(x_rung,x_func)

!------------------------------------------!
! CORRELATION: read rung of Jacob's ladder !
!------------------------------------------!

  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) c_rung,c_func

  select case (c_rung) ! correlation functionals
 
    case (0) ! Hartree

      select case (c_func)

        case ('H')
     
          c_DFA = 1
     
        case default
     
          call print_warning('!!! Hartree correlation functional not available !!!')
          stop

      end select

    case (1) ! LDA
   
      select case (c_func)
 
        case ('W38')
  
          c_DFA = 1
  
        case ('PW92')
  
          c_DFA = 2
  
        case ('VWN3')
  
          c_DFA = 3
  
        case ('VWN5')
  
          c_DFA = 4
  
        case ('eVWN5')
  
          c_DFA = 5
  
        case default
 
          call print_warning('!!! LDA correlation functional not available !!!')
          stop
 
      end select

    case (2) ! GGA
   
    select case (c_func)

      case ('LYP')
 
        c_DFA = 1
 
      case ('PBE')
 
        c_DFA = 2
 
      case default

        call print_warning('!!! GGA correlation functional not available !!!')
        stop

    end select

    case (3) ! MGGA
   
      select case (c_func)
 
        case default
 
          call print_warning('!!! MGGA correlation functional not available !!!')
          stop
 
      end select

    case (4) ! Hybrid
   
      select case (c_func)
 
        case ('HF')
  
          c_DFA = 1
  
        case ('B3LYP')
  
          c_DFA = 2
  
        case ('BHHLYP')
 
          c_DFA = 3
  
        case ('PBE')
  
          c_DFA = 4
  
        case default
 
          call print_warning('!!! Hybrid correlation functional not available !!!')
          stop

      end select

    case default 

      call print_warning('!!! Correlation rung not available !!!')
      stop

  end select

! Select rung for correlation

  write(*,*)
  write(*,*) '*******************************************************************'
  write(*,*) '*                       Correlation rung                          *'
  write(*,*) '*******************************************************************'

  call select_rung(c_rung,c_func)

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

  do iEns=1,maxEns
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

  allocate(nEl(maxEns))
  nEl(:) = 0d0
  do iEns=1,maxEns
    do iBas=1,nBas
      nEl(iEns) = nEl(iEns) + occnum(iBas,1,iEns) + occnum(iBas,2,iEns)
    end do
  end do

  doNcentered = .false.

  read(1,*)
  read(1,*) (wEns(iEns),iEns=2,nEns)
  read(1,*)
  read(1,*) answer

  if(answer == 'T') doNcentered = .true.

  wEns(1) = 1d0
  do iEns=2,nEns
    wEns(1) = wEns(1) - wEns(iEns)
  end do

  if (doNcentered) then

    do iEns=2,nEns
      wEns(iEns) = (nEl(1)/nEl(iEns))*wEns(iEns)
    end do

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
  write(*,*)' parameters for w1-dependent exchange functional coefficient '
  write(*,*)'----------------------------------------------------------'
  call matout(3,1,aCC_w1)
  write(*,*)

  write(*,*)'----------------------------------------------------------'
  write(*,*)' parameters for w2-dependent exchange functional coefficient '
  write(*,*)'----------------------------------------------------------'
  call matout(3,1,aCC_w2)
  write(*,*) 

! Close file with options

  close(unit=1)

end subroutine read_options_dft
