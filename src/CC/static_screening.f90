subroutine static_screening(nBas,nC,nO,nV,nR,eW,ERI,dbERI)

! Compute the four-index tensor of the static screening W 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR

  double precision,intent(in)    :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)    :: eW(nBas)

! Local variables

  double precision              :: EcRPA
  double precision              :: eta
  logical                       :: TDA
  double precision              :: chi
  integer                       :: ispin
  integer                       :: nS
  integer                       :: p,q,r,s
  integer                       :: m

  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

! Output variables

  double precision,intent(out)   :: dbERI(nBas,nBas,nBas,nBas)

! Initialize 

  nS = (nO - nC)*(nV - nR)

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  ispin = 3
  EcRPA = 0d0

  eta = 0d0
  TDA = .false.

  call phLR(ispin,.true.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,EcRPA,Om,XpY,XmY)
  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,rho)

  do p=1,nBas
    do q=1,nBas
      do r=1,nBas
        do s=1,nBas
          chi = 0d0
          do m=1,nS
            chi = chi + rho(p,s,m)*rho(q,r,m)/Om(m)
          enddo

          dbERI(p,q,r,s)= ERI(p,q,r,s) - ERI(p,q,s,r) + 2d0*chi

        enddo
      enddo
    enddo
  enddo

end subroutine 
