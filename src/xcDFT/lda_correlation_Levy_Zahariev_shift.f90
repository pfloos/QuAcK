subroutine lda_correlation_Levy_Zahariev_shift(DFA,nEns,wEns,nGrid,weight,rho,EcLZ)

! Compute the lda correlation part of the Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Output variables

  double precision,intent(out)  :: EcLZ(nsp)

! Select correlation functional

  select case (DFA)

!   Wigner's LDA correlation functional: Wigner, Trans. Faraday Soc. 34 (1938) 678

    case ('W38')

      call W38_lda_correlation_Levy_Zahariev_shift(nGrid,weight(:),rho(:,:),EcLZ(:))

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('VWN5')

      call VWN5_lda_correlation_Levy_Zahariev_shift(nGrid,weight(:),rho(:,:),EcLZ(:))

!   Loos-Fromager weight-dependent correlation functional: Loos and Fromager (in preparation)

    case ('LF19')

      call LF19_lda_correlation_Levy_Zahariev_shift(nEns,wEns,nGrid,weight(:),rho(:,:),EcLZ(:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select
 
end subroutine lda_correlation_Levy_Zahariev_shift
