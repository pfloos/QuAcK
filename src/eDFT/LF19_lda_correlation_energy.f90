subroutine LF19_lda_correlation_energy(nEns,wEns,nGrid,weight,rho,Ec)

! Compute eLDA correlation energy 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  logical                       :: LDA_centered = .false.
  integer                       :: iEns,isp
  double precision              :: EcLDA(nsp)
  double precision,allocatable  :: aLF(:,:)
  double precision,allocatable  :: EceLDA(:,:)

! Output variables

  double precision              :: Ec(nsp)

! Allocation

  allocate(aLF(3,nEns),EceLDA(nsp,nEns))

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

    call elda_correlation_energy(nEns,aLF(:,iEns),nGrid,weight(:),rho(:,:),EceLDA(:,iEns))

  end do

! LDA-centered functional

  EcLDA(:) = 0d0

  if(LDA_centered) then 

    call VWN5_lda_correlation_energy(nGrid,weight(:),rho(:,:),EcLDA(:))

    do iEns=1,nEns
      do isp=1,nsp

        EceLDA(isp,iEns) = EceLDA(isp,iEns) + EcLDA(isp) - EceLDA(isp,1)

      end do
    end do

  end if

! Weight-denpendent functional for ensembles

  Ec(:) = 0d0

  do iEns=1,nEns
    do isp=1,nsp

      Ec(isp) = Ec(isp) + wEns(iEns)*EceLDA(isp,iEns)

    end do
  end do

end subroutine LF19_lda_correlation_energy
