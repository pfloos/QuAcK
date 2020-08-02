subroutine unrestricted_lda_exchange_individual_energy(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,rho,Ex,Cx_choice)

! Compute LDA exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice

! Output variables

  double precision              :: Ex

! Select correlation functional

  select case (DFA)

    case ('S51')

      call US51_lda_exchange_individual_energy(nGrid,weight,rhow,rho,Ex)

    case ('CC')

      call UCC_lda_exchange_individual_energy(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,rho,Ex,Cx_choice)

    case default

      call print_warning('!!! LDA exchange individual energy not available !!!')
      stop

  end select

end subroutine unrestricted_lda_exchange_individual_energy
