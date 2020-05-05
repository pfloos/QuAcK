subroutine lda_exchange_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,ExDD)

! Compute the exchange LDA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Select correlation functional

  select case (DFA)

    case ('S51')

      ExDD(:) = 0d0

    case ('RS51')

      ExDD(:) = 0d0

    case ('RMFL20')

      call RMFL20_lda_exchange_derivative_discontinuity(nEns,wEns,nGrid,weight(:),rhow(:),ExDD(:))

    case ('RCC')

      call RCC_lda_exchange_derivative_discontinuity(nEns,wEns,nGrid,weight(:),rhow(:),ExDD(:))

    case default

      call print_warning('!!! LDA exchange derivative discontinuity not available !!!')
      stop

  end select
 
end subroutine lda_exchange_derivative_discontinuity
