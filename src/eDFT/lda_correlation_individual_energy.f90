subroutine lda_correlation_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,rho,Ec)

! Compute LDA correlation energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Output variables

  double precision              :: Ec(nsp)

! Select correlation functional

  select case (DFA)

!   Wigner's LDA correlation functional: Wigner, Trans. Faraday Soc. 34 (1938) 678

    case ('W38')

      call W38_lda_correlation_individual_energy(nGrid,weight(:),rhow(:,:),rho(:,:),Ec(:))

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('VWN5')

      call VWN5_lda_correlation_individual_energy(nGrid,weight(:),rhow(:,:),rho(:,:),Ec(:))

!   Loos-Fromager weight-dependent correlation functional: Loos and Fromager (in preparation)

    case ('MFL20')

      call MFL20_lda_correlation_individual_energy(nEns,wEns,nGrid,weight(:),rhow(:,:),rho(:,:),Ec(:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine lda_correlation_individual_energy
