subroutine exchange_energy(rung,DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P,FxHF, &
                                        rho,drho,Cx_choice,doNcentered,Ex)

! Compute the exchange energy

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
  double precision,intent(in)   :: FxHF(nBas,nBas)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Local variables

! Output variables

  double precision,intent(out)  :: Ex

  select case (rung)

!   Hartree calculation

    case(0) 

      Ex = 0d0

!   LDA functionals

    case(1) 

      call lda_exchange_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,rho,Cx_choice,doNcentered,Ex)

!   GGA functionals

    case(2) 

      call gga_exchange_energy(DFA,nEns,wEns,nCC,aCC,nGrid,weight,rho,drho,Cx_choice,Ex)

!   MGGA functionals

    case(3) 

      call mgga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

!   Hybrid functionals

    case(4) 

      call hybrid_exchange_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P,FxHF, &
                                        rho,drho,Cx_choice,doNcentered,Ex)

  end select
 
end subroutine exchange_energy
