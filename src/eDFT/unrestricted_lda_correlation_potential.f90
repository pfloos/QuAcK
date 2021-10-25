subroutine unrestricted_lda_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,rho,Fc)

! Select LDA correlation potential

  implicit none
include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
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

    case (1)

      call UW38_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case (2)

      call UPW92_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case (3)

      call UVWN3_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case (4)

      call UVWN5_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case (5)

      call UC16_lda_correlation_potential(nGrid,weight(:),nBas,AO(:,:),rho(:,:),Fc(:,:,:))

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine unrestricted_lda_correlation_potential
