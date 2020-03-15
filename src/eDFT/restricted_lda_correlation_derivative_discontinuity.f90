subroutine restricted_lda_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

! Compute the correlation LDA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables

  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec(nEns)

! Select correlation functional

  select case (DFA)

!   Wigner's LDA correlation functional: Wigner, Trans. Faraday Soc. 34 (1938) 678

    case ('W38')

      Ec(:) = 0d0

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('VWN5')

      Ec(:) = 0d0

!   Loos-Fromager weight-dependent correlation functional: Loos and Fromager (in preparation)

    case ('LF19')

      call RMFL20_lda_correlation_derivative_discontinuity(nEns,wEns,nGrid,weight(:),rhow(:),Ec(:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select
 
end subroutine restricted_lda_correlation_derivative_discontinuity
