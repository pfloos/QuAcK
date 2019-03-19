subroutine LF19_lda_correlation_potential(nEns,wEns,nGrid,weight,nBas,AO,rho,Fc)

! Compute Loos-Fromager weight-dependent LDA correlation potential

  implicit none
include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  logical                       :: LDA_centered = .false.
  integer                       :: iEns
  double precision,allocatable  :: aLF(:,:)
  double precision,allocatable  :: FcLDA(:,:,:)
  double precision,allocatable  :: FceLDA(:,:,:,:)

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Allocation

  allocate(aLF(3,nEns),FcLDA(nBas,nBas,nspin),FceLDA(nBas,nBas,nspin,nEns))

! Parameters for weight-dependent LDA correlation functional

  aLF(1,1) = -0.0238184d0
  aLF(2,1) = +0.00575719d0
  aLF(3,1) = +0.0830576d0

  aLF(1,2) = -0.0282814d0
  aLF(2,2) = +0.00340758d0
  aLF(3,2) = +0.0663967d0

  aLF(1,3) = -0.0144633d0
  aLF(2,3) = -0.0504501d0
  aLF(3,3) = +0.0331287d0

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

    call elda_correlation_potential(nEns,aLF(:,iEns),nGrid,weight,nBas,AO,rho,FceLDA(:,:,:,iEns))

  end do

! LDA-centered functional

  FcLDA(:,:,:) = 0d0

  if(LDA_centered) then 

    call VWN5_lda_correlation_potential(nGrid,weight,nBas,AO,rho,FcLDA)

    do iEns=1,nEns

      FceLDA(:,:,:,iEns) = FceLDA(:,:,:,iEns) + FcLDA(:,:,:) - FceLDA(:,:,:,1)

    end do

  end if

! Weight-denpendent functional for ensembles

  Fc(:,:,:) = 0d0

  do iEns=1,nEns

    Fc(:,:,:) = Fc(:,:,:) + wEns(iEns)*FceLDA(:,:,:,iEns)

  enddo

end subroutine LF19_lda_correlation_potential
