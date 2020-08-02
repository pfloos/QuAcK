subroutine restricted_lda_exchange_potential(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,AO,rho,Fx,Cx_choice)

! Select LDA correlation potential

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
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Select exchange functional

  select case (DFA)

    case ('S51')

      call RS51_lda_exchange_potential(nGrid,weight,nBas,AO,rho,Fx)

    case ('MFL20')

      call RMFL20_lda_exchange_potential(LDA_centered,nEns,wEns,nGrid,weight,nBas,AO,rho,Fx)

    case ('CC')

      call RCC_lda_exchange_potential(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,AO,rho,Fx,Cx_choice)

    case default

      call print_warning('!!! LDA exchange functional not available !!!')
      stop

  end select

end subroutine restricted_lda_exchange_potential
