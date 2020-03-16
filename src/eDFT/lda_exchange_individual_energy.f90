subroutine lda_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute LDA exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision              :: Ex

! Select correlation functional

  select case (DFA)

    case ('S51')

!     call S51_lda_exchange_individual_energy(nGrid,weight(:),rhow(:),rho(:),Ex)

    case ('MFL20')

!     call MFL20_lda_exchange_individual_energy(nEns,wEns,nGrid,weight(:),rhow(:),rho(:),Ex)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine lda_exchange_individual_energy
