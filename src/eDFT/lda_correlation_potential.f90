subroutine lda_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,rho,Fc)

! Select LDA correlation potential

  implicit none
include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Select correlation functional

  select case (DFA)

!   Hartree-Fock 

    case ('HF')

      Fc(:,:,:) = 0d0

!   Wigner's LDA correlation functional: Wigner, Trans. Faraday Soc. 34 (1938) 678

    case ('W38')

      call W38_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

!   Vosko, Wilk and Nusair's functional V: Can. J. Phys. 58 (1980) 1200

    case ('VWN5')

      call VWN5_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

!   Chachiyo's LDA correlation functional: Chachiyo, JCP 145 (2016) 021101

    case ('C16')

      call C16_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

!   Restricted version of the Marut-Fromager-Loos weight-dependent correlation functional

    case ('RMFL20')

      call RMFL20_lda_correlation_potential(nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine lda_correlation_potential
