subroutine RMFL20_lda_correlation_individual_energy(LDA_centered,nEns,wEns,nGrid,weight,rhow,rho,Ec)

! Compute eLDA correlation energy 

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
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

  aMFL(1,2) = -0.0282814d0
  aMFL(2,2) = +0.00273925d0
  aMFL(3,2) = +0.0664914d0

  aMFL(1,3) = -0.0144633d0
  aMFL(2,3) = -0.0506019d0
  aMFL(3,3) = +0.0331417d0

! Compute correlation energy for ground- and doubly-excited states

  do iEns=1,nEns
    call restricted_elda_correlation_individual_energy(aMFL(1:3,iEns),nGrid,weight(:),rhow(:),rho(:),EceLDA(iEns))
  end do

! LDA-centered functional

  if(LDA_centered) then 
  
    call RVWN5_lda_correlation_individual_energy(nGrid,weight(:),rhow(:),rho(:),EcLDA)

    EceLDA(3) = EcLDA + wEns(3)*(EceLDA(3) - EceLDA(1))
    EceLDA(2) = EcLDA + wEns(2)*(EceLDA(2) - EceLDA(1))
    EceLDA(1) = EcLDA

  end if

! Weight-denpendent functional for ensembles

  Ec = dot_product(wEns(:),EceLDA(:))

end subroutine RMFL20_lda_correlation_individual_energy
