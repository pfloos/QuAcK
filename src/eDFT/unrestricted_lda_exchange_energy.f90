subroutine unrestricted_lda_exchange_energy(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rho,Ex)

! Select LDA exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision,intent(out)  :: Ex

! Select correlation functional

  select case (DFA)

    case ('S51')

      call US51_lda_exchange_energy(nGrid,weight,rho,Ex)

    case ('CC')

      call UCC_lda_exchange_energy(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rho,Ex)

    case default

      call print_warning('!!! LDA exchange functional not available !!!')
      stop

  end select

end subroutine unrestricted_lda_exchange_energy
