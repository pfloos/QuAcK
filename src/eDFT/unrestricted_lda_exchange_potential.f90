subroutine unrestricted_lda_exchange_potential(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,rho, &
                                               Cx_choice,doNcentered,Fx)

! Select LDA correlation potential

  implicit none

  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Select exchange functional

  select case (DFA)

    case (1)

      call US51_lda_exchange_potential(nGrid,weight,nBas,AO,rho,Fx)

    case (2)

      call UCC_lda_exchange_potential(nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,rho,Cx_choice,doNcentered,Fx)

    case default

      call print_warning('!!! LDA exchange potential not available !!!')
      stop

  end select

end subroutine unrestricted_lda_exchange_potential
