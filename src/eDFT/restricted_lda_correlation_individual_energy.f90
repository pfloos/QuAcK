subroutine restricted_lda_correlation_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,rho,Ec)

! Compute LDA correlation energy for individual states

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

  double precision              :: Ec

! Select correlation functional

  select case (DFA)

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('RVWN5')

      call RVWN5_lda_correlation_individual_energy(nGrid,weight(:),rhow(:),rho(:),Ec)

!   marut-Fromager-Loos weight-dependent correlation functional

    case ('RMFL20')

      call RMFL20_lda_correlation_individual_energy(nEns,wEns,nGrid,weight(:),rhow(:),rho(:),Ec)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine restricted_lda_correlation_individual_energy
