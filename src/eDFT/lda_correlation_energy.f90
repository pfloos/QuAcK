subroutine lda_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,Ec)

! Select LDA correlation functional

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

  double precision,intent(out)  :: Ec(nsp)

! Select correlation functional

  select case (DFA)

!   Hartree-Fock

    case ('HF')

      Ec(:) = 0d0

    case ('W38')

      call W38_lda_correlation_energy(nGrid,weight(:),rho(:,:),Ec(:))

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('VWN5')

      call VWN5_lda_correlation_energy(nGrid,weight(:),rho(:,:),Ec(:))

!   Chachiyo's LDA correlation functional: Chachiyo, JCP 145 (2016) 021101

    case ('C16')

      call C16_lda_correlation_energy(nGrid,weight(:),rho(:,:),Ec(:))

!   Restricted version of the Marut-Fromager-Loos weight-dependent correlation functional

    case ('RMFL20')

      call RMFL20_lda_correlation_energy(nEns,wEns(:),nGrid,weight(:),rho(:,:),Ec(:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine lda_correlation_energy
