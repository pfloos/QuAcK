subroutine lda_exchange_energy(DFA,LDA_centered,nEns,wEns,nGrid,weight,rho,Ex)

! Select LDA exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision,intent(out)  :: Ex

! Select correlation functional

  select case (DFA)

    case ('S51')

      call S51_lda_exchange_energy(nGrid,weight,rho,Ex)

    case ('RS51')

      call RS51_lda_exchange_energy(nGrid,weight,rho,Ex)

    case ('RMFL20')

      call RMFL20_lda_exchange_energy(LDA_centered,nEns,wEns,nGrid,weight,rho,Ex)

    case ('RCC')

      call RCC_lda_exchange_energy(nEns,wEns,nGrid,weight,rho,Ex)

    case default

      call print_warning('!!! LDA exchange functional not available !!!')
      stop

  end select

end subroutine lda_exchange_energy
