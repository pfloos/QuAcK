subroutine unrestricted_exchange_potential(rung,DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,P, & 
                                           ERI,AO,dAO,rho,drho,Fx,FxHF,Cx_choice,doNcentered)

! Compute the exchange potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: DFA
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
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

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

      call unrestricted_lda_exchange_potential(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,AO,rho,Fx,&
                                               Cx_choice,doNcentered)

!   GGA functionals

    case(2) 

      call unrestricted_gga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_exchange_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_exchange_potential(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,P, &
                                                  ERI,AO,dAO,rho,drho,Fx,FxHF,Cx_choice)

  end select

end subroutine unrestricted_exchange_potential
