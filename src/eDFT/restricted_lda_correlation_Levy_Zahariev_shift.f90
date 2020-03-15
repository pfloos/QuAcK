subroutine restricted_lda_correlation_Levy_Zahariev_shift(DFA,nEns,wEns,nGrid,weight,rho,EcLZ)

! Compute the lda correlation part of the Levy-Zahariev shift

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

  double precision,intent(out)  :: EcLZ

! Select correlation functional

  select case (DFA)

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('RVWN5')

      call RVWN5_lda_correlation_Levy_Zahariev_shift(nGrid,weight(:),rho(:),EcLZ)

!   Marut-Fromager-Loos weight-dependent correlation functional

    case ('RMFL20')

      call RMFL20_lda_correlation_Levy_Zahariev_shift(nEns,wEns,nGrid,weight(:),rho(:),EcLZ)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select
 
end subroutine restricted_lda_correlation_Levy_Zahariev_shift
