subroutine mgga_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

! Compute the correlation MGGA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)

! Local variables

! Output variables

  double precision,intent(out)  :: Ec(nsp,nEns)

! Select correlation functional

  select case (DFA)

    case default

      call print_warning('!!! MGGA correlation functional not available !!!')
      stop

  end select
 
end subroutine mgga_correlation_derivative_discontinuity
