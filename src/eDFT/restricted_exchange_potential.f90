subroutine restricted_exchange_potential(rung,DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,P, & 
                                         ERI,AO,dAO,rho,drho,Fx,FxHF)

! Compute the exchange potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision,allocatable  :: FxLDA(:,:),FxGGA(:,:)
  double precision              :: cX,aX

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas),FxHF(nBas,nBas)

! Memory allocation

  select case (rung)

!   Hartree calculation

    case(0) 

      Fx(:,:) = 0d0

!   LDA functionals

    case(1) 

      call restricted_lda_exchange_potential(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,AO,rho,Fx)

!   GGA functionals

    case(2) 

      call restricted_gga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

!   Hybrid functionals

    case(4) 

      allocate(FxLDA(nBas,nBas),FxGGA(nBas,nBas))

      cX = 0.20d0
      aX = 0.72d0

      call restricted_lda_exchange_potential(DFA,nGrid,weight,nBas,AO,rho,FxLDA)
      call restricted_gga_exchange_potential(DFA,nGrid,weight,nBas,AO,dAO,rho,drho,FxGGA)
      call restricted_fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = FxLDA(:,:)                   &
              + cX*(FxHF(:,:)  - FxLDA(:,:)) &
              + aX*(FxGGA(:,:) - FxLDA(:,:))  

!   Hartree-Fock calculation

    case(666) 

      call restricted_fock_exchange_potential(nBas,P,ERI,FxHF)

      Fx(:,:) = FxHF(:,:)

  end select

end subroutine restricted_exchange_potential
