subroutine lda_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

! Compute the correlation LDA part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)

! Local variables

  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec(nsp,nEns)

! Select correlation functional

  select case (DFA)

    case (1)

      Ec(:,:) = 0d0

    case (2)

      Ec(:,:) = 0d0

    case (3)

      Ec(:,:) = 0d0

    case (4)

      Ec(:,:) = 0d0

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select
 
end subroutine lda_correlation_derivative_discontinuity
