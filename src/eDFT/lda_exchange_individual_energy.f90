subroutine lda_exchange_individual_energy(DFA,LDA_centered,nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute LDA exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
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

    case ('RS51')

      call RS51_lda_exchange_individual_energy(nGrid,weight(:),rhow(:),rho(:),Ex)

    case ('RMFL20')

      call RMFL20_lda_exchange_individual_energy(LDA_centered,nEns,wEns,nGrid,weight(:),rhow(:),rho(:),Ex)

    case ('RGIC')

      call RGIC_lda_exchange_individual_energy(nEns,wEns,nGrid,weight(:),rhow(:),rho(:),Ex)

    case default

      call print_warning('!!! LDA exchange individual energy not available !!!')
      stop

  end select

end subroutine lda_exchange_individual_energy
