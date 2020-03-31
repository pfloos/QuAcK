subroutine restricted_lda_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,rho,Fc)

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
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas)

! Select correlation functional

  select case (DFA)

!   Hartree-Fock 

    case ('HF')

      Fc(:,:) = 0d0

    case ('RVWN5')

      call RVWN5_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:),Fc(:,:))

    case ('RMFL20')

      call RMFL20_lda_correlation_potential(nEns,wEns(:),nGrid,weight(:),nBas,AO(:,:),rho(:),Fc(:,:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine restricted_lda_correlation_potential
