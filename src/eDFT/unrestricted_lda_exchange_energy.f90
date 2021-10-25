subroutine unrestricted_lda_exchange_energy(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rho,Ex,Cx_choice)

! Select LDA exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice

! Output variables

  double precision,intent(out)  :: Ex

! Select correlation functional

  select case (DFA)

    case (1)

      call US51_lda_exchange_energy(nGrid,weight,rho,Ex)

    case (2)

      call UCC_lda_exchange_energy(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rho,Ex,Cx_choice)

    case default

      call print_warning('!!! LDA exchange energy not available !!!')
      
      stop

  end select

end subroutine unrestricted_lda_exchange_energy
