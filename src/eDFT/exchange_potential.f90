subroutine exchange_potential(rung,DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P, & 
                                           ERI,AO,dAO,rho,drho,Cx_choice,doNcentered,Fx,FxHF)

! Compute the exchange potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
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
  double precision              :: cX,aX

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)
  double precision,intent(out)  :: FxHF(nBas,nBas)

! Memory allocation

  select case (rung)

!   Hartree calculation

    case(0) 

      Fx(:,:) = 0d0

!   LDA functionals

    case(1) 

      call lda_exchange_potential(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,rho,&
                                               Cx_choice,doNcentered,Fx)

!   GGA functionals

    case(2) 

      call gga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

!   MGGA functionals

    case(3) 

      call mgga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

!   Hybrid functionals

    case(4) 

      call hybrid_exchange_potential(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P, &
                                                  ERI,AO,dAO,rho,drho,Cx_choice,doNcentered,Fx,FxHF)

  end select

end subroutine exchange_potential
