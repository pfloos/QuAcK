subroutine unrestricted_mgga_exchange_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,drhow,ExDD)

! Compute the exchange MGGA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)

! Local variables


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Select exchange functional

  select case (DFA)

    case default

      call print_warning('!!! MGGA exchange derivative discontinuity not available !!!')
      stop

  end select
 
end subroutine unrestricted_mgga_exchange_derivative_discontinuity
