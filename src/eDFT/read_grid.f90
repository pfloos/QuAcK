subroutine read_grid(SGn,nRad,nAng,nGrid)

! Read grid type

  implicit none

! Input variables

  integer,intent(in)            :: SGn

! Output variables

  integer,intent(out)           :: nRad
  integer,intent(out)           :: nAng
  integer,intent(out)           :: nGrid

  write(*,*)'----------------------------------------------------------'
  write(*,'(A22,I1)')'  Quadrature grid: SG-',SGn
  write(*,*)'----------------------------------------------------------'

  select case (SGn)

    case(0)
     nRad = 23
     nAng = 170

    case(1)
      nRad = 50
      nAng = 194

    case(2)
      nRad = 75
      nAng = 302

    case(3)
      nRad = 99
      nAng = 590

    case default
      call print_warning('!!! Quadrature grid not available !!!')
      stop

  end select

  nGrid = nRad*nAng

end subroutine read_grid
