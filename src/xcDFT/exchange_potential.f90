subroutine exchange_potential(rung,nGrid,weight,nBas,P,ERI,AO,dAO,rho,drho,Fx,FxHF)

! Compute the exchange potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  double precision,allocatable  :: FxLDA(:,:),FxGGA(:,:)
  double precision              :: cX,aX,aC

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas),FxHF(nBas,nBas)

! Memory allocation

  allocate(FxLDA(nBas,nBas),FxGGA(nBas,nBas))

  FxLDA(:,:) = 0d0
  FxGGA(:,:) = 0d0

  select case (rung)

!   Hartree calculation
    case(0) 

      Fx(:,:) = 0d0

!   LDA functionals
    case(1) 

      call lda_exchange_potential(nGrid,weight,nBas,AO,rho,FxLDA)

      Fx(:,:) = FxLDA(:,:)

!   GGA functionals
    case(2) 

      call gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)

      Fx(:,:) = FxGGA(:,:)

!   Hybrid functionals
    case(4) 

      cX = 0.20d0
      aX = 0.72d0
      aC = 0.81d0

      call lda_exchange_potential(nGrid,weight,nBas,AO,rho,FxLDA)
      call gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)
      call fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = FxLDA(:,:)                   &
              + cX*(FxHF(:,:)  - FxLDA(:,:)) &
              + aX*(FxGGA(:,:) - FxLDA(:,:))  

!   Hartree-Fock calculation
    case(666) 

      call fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = FxHF(:,:)

  end select

end subroutine exchange_potential
