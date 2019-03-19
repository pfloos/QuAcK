subroutine LF19_lda_correlation_derivative_discontinuity(nEns,wEns,nGrid,weight,rhow,Ec)

! Compute eLDA correlation part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)

! Local variables

  integer                       :: iEns,jEns
  double precision,allocatable  :: aLF(:,:)
  double precision              :: dEc(nsp,nEns)
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: Ec(nsp,nEns)

! Allocation

  allocate(aLF(3,nEns))

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

    call elda_correlation_energy(nEns,aLF(:,iEns),nGrid,weight(:),rhow(:,:),dEc(:,iEns))

  end do

  Ec(:,:) = 0d0

  do iEns=1,nEns
    do jEns=1,nEns

      Ec(:,iEns) = Ec(:,iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*(dEc(:,jEns) - dEc(:,1))

    end do
  end do

end subroutine LF19_lda_correlation_derivative_discontinuity
