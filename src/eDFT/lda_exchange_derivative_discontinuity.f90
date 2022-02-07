subroutine lda_exchange_derivative_discontinuity(DFA,nEns,wEns,nCC,aCC,nGrid,weight,rhow,&
                                                        Cx_choice,doNcentered,kappa,ExDD)

! Compute the exchange LDA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
 
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered
  double precision,intent(in)   :: kappa(nEns)

! Local variables


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Select exchange functional

  select case (DFA)

    case (1)

      ExDD(:) = 0d0

    case (2)

      call CC_lda_exchange_derivative_discontinuity(nEns,wEns,nCC,aCC,nGrid,weight,rhow,&
                                                     Cx_choice,doNcentered,kappa,ExDD)

    case default

      call print_warning('!!! LDA exchange derivative discontinuity not available !!!')
      stop

  end select
 
end subroutine lda_exchange_derivative_discontinuity
