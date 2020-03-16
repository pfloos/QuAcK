subroutine MFL20_lda_exchange_Levy_Zahariev_shift(nEns,wEns,nGrid,weight,rho,ExLZ)

! Compute the Marut-Fromager-Loos LDA exchange contribution to Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  logical                       :: LDA_centered = .true.
  integer                       :: iEns
  double precision              :: ExLZLDA
  double precision,allocatable  :: aMFL(:,:)
  double precision,allocatable  :: ExLZeLDA(:)

! Output variables

  double precision,intent(out)  :: ExLZ


! Allocation

  allocate(aMFL(3,nEns),ExLZeLDA(nEns))

! Parameters for weight-dependent LDA correlation functional

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

!   call elda_exchange_Levy_Zahariev_shift(nEns,aMFL(:,iEns),nGrid,weight(:),rho(:),ExLZeLDA(iEns))

  end do

! LDA-centered functional

  ExLZLDA = 0d0

  if(LDA_centered) then

!   call S51_lda_exchange_Levy_Zahariev_shift(nGrid,weight(:),rho(:),ExLZLDA)

    do iEns=1,nEns
 
      ExLZeLDA(iEns) = ExLZeLDA(iEns) + ExLZLDA - ExLZeLDA(1)
 
    end do

  end if

! Weight-denpendent functional for ensembles

  ExLZ = 0d0

  do iEns=1,nEns

    ExLZ = ExLZ + wEns(iEns)*ExLZeLDA(iEns)

  enddo

end subroutine MFL20_lda_exchange_Levy_Zahariev_shift
