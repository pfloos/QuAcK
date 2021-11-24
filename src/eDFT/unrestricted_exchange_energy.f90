subroutine unrestricted_exchange_energy(rung,DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P,FxHF, &
                                        rho,drho,Ex,Cx_choice,doNcentered)

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

      call unrestricted_lda_exchange_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,&
                                            rho,Ex,Cx_choice,doNcentered)

!   GGA functionals

    case(2) 

      call unrestricted_gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ex)

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_exchange_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,P,FxHF, &
                                        rho,drho,Ex,Cx_choice)

  end select
 
end subroutine unrestricted_exchange_energy
