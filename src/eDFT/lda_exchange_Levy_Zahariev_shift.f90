subroutine lda_exchange_Levy_Zahariev_shift(DFA,nEns,wEns,nGrid,weight,rho,ExLZ)

! Compute the exchange LDA part of the Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision,intent(out)  :: ExLZ

! Select correlation functional

  select case (DFA)

    case ('S51')

      call print_warning('!!! Exchange Levzy-Zahariev shift NYI for S51 !!!')
      stop

    case ('MFL20')

      call MFL20_lda_exchange_Levy_Zahariev_shift(nEns,wEns,nGrid,weight(:),rho(:),ExLZ)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select
 
end subroutine lda_exchange_Levy_Zahariev_shift
