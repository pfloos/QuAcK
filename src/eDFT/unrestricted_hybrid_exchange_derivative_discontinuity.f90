subroutine unrestricted_hybrid_exchange_derivative_discontinuity(DFA,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,&
                                                              Cx_choice,doNcentered,kappa,ExDD)

! Compute the exchange part of the derivative discontinuity for hybrid functionals

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
 
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

    case ('HF')

      ExDD(:) = 0d0

    case ('B3')

      ExDD(:) = 0d0

    case ('PBE')

      ExDD(:) = 0d0

    case default

      call print_warning('!!! Hybrid exchange derivative discontinuity not available !!!')
      stop

  end select
 
end subroutine unrestricted_hybrid_exchange_derivative_discontinuity
