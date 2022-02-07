subroutine hybrid_exchange_potential(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P, & 
                                                  ERI,AO,dAO,rho,drho,Cx_choice,doNcentered,Fx,FxHF)

! Compute the exchange potential for hybrid functionals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Local variables

  double precision,allocatable  :: FxLDA(:,:)
  double precision,allocatable  :: FxGGA(:,:)
  double precision              :: a0
  double precision              :: aX

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)
  double precision,intent(out)  :: FxHF(nBas,nBas)

! Memory allocation

  select case (DFA)

    case(1) 

      call fock_exchange_potential(nBas,P,ERI,FxHF)
      Fx(:,:) = FxHF(:,:)

    case(2) 

      allocate(FxLDA(nBas,nBas),FxGGA(nBas,nBas))

      a0 = 0.20d0
      aX = 0.72d0

      call lda_exchange_potential(1,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight, & 
                                               nBas,AO,rho,Cx_choice,doNcentered,FxLDA)
      call gga_exchange_potential(2,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)
      call fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = FxLDA(:,:)                   &
              + a0*(FxHF(:,:)  - FxLDA(:,:)) &
              + aX*(FxGGA(:,:) - FxLDA(:,:))  

    case(3) 

      allocate(FxGGA(nBas,nBas))

      call gga_exchange_potential(2,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)
      call fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = 0.5d0*FxHF(:,:) + 0.5d0*FxGGA(:,:)

    case(4) 

      allocate(FxGGA(nBas,nBas))

      call gga_exchange_potential(3,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)
      call fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = 0.25d0*FxHF(:,:) + 0.75d0*FxGGA(:,:)

    case default

      call print_warning('!!! Hybrid exchange potential not available !!!')
      stop

  end select

end subroutine hybrid_exchange_potential
