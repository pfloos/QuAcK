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

!   Slater's LDA correlation functional: Slater, Phys. Rev. 81 (1951) 385

    case ('S51')

      call S51_lda_exchange_energy(nGrid,weight,rho,Ex)

    case default

      call print_warning('!!! LDA exchange functional not available !!!')
      stop

  end select

end subroutine lda_exchange_energy
