subroutine read_grid(SGn,radial_precision,nRad,nAng)

! Read grid type

  implicit none

! Input variables

  integer,intent(in)            :: SGn

! Output variables

  double precision,intent(out)  :: radial_precision
  integer,intent(out)           :: nRad
  integer,intent(out)           :: nAng

  write(*,*)'----------------------------------------------------------'
  write(*,'(A22,I1)')'  Quadrature grid: SG-',SGn
  write(*,*)'----------------------------------------------------------'

  select case (SGn)

    case(0)
      radial_precision = 1d-5
      nRad = 23
      nAng = 170

    case(1)
      radial_precision = 1d-7
      nRad = 50
      nAng = 194

    case(2)
      radial_precision = 1d-9
      nRad = 75
      nAng = 302

    case(3)
      radial_precision = 1d-11
      nRad = 99
      nAng = 590

    case default
      call print_warning('!!! Quadrature grid not available !!!')
      stop

  end select

end subroutine read_grid
