subroutine lda_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,Ex)

! Select LDA exchange functional

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

  double precision,intent(out)  :: Ex

! Select correlation functional

  select case (DFA)

!   Slater's LDA correlation functional

    case ('RS51')

      call RS51_lda_exchange_energy(nGrid,weight,rho,Ex)

!   Restricted version of Slater's LDA correlation functional

    case ('S51')

      call S51_lda_exchange_energy(nGrid,weight,rho,Ex)

!   Restricted version of the weight-dependent Marut-Fromager-Loos 2020 exchange functional

    case ('RMFL20')

      call RMFL20_lda_exchange_energy(nEns,wEns,nGrid,weight,rho,Ex)

    case default

      call print_warning('!!! LDA exchange functional not available !!!')
      stop

  end select

end subroutine lda_exchange_energy
