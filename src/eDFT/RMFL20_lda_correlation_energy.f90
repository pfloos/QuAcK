subroutine RMFL20_lda_correlation_energy(nEns,wEns,nGrid,weight,rho,Ec)

! Compute the restricted version of the Marut-Fromager-Loos weight-dependent correlation functional
! The RMFL20 is a two-state, single-weight correlation functional for spin-unpolarized systems

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iEns
  double precision              :: EcLDA
  double precision,allocatable  :: aMFL(:,:)
  double precision,allocatable  :: EceLDA(:)

! Output variables

  double precision              :: Ec

! Allocation

  allocate(aMFL(3,nEns),EceLDA(nEns))

! Parameters for weight-dependent LDA correlation functional

  aMFL(1,1) = -0.0238184d0
  aMFL(2,1) = +0.00540994d0
  aMFL(3,1) = +0.0830766d0

  aMFL(1,2) = -0.0144633d0
  aMFL(2,2) = -0.0506019d0
  aMFL(3,2) = +0.0331417d0

! Compute correlation energy for ground and doubly-excited states

  do iEns=1,nEns

    call restricted_elda_correlation_energy(nEns,aMFL(:,iEns),nGrid,weight(:),rho(:),EceLDA(iEns))

  end do

! LDA-centered functional

  EcLDA = 0d0

  call RVWN5_lda_correlation_energy(nGrid,weight(:),rho(:),EcLDA)

  do iEns=1,nEns

    EceLDA(iEns) = EceLDA(iEns) + EcLDA - EceLDA(1)

  end do

! Weight-denpendent functional for ensembles

  Ec = 0d0

  do iEns=1,nEns

    Ec = Ec + wEns(iEns)*EceLDA(iEns)

  end do

end subroutine RMFL20_lda_correlation_energy
