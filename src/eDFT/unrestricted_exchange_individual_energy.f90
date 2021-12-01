subroutine unrestricted_exchange_individual_energy(rung,DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas, & 
                                                   ERI,Pw,rhow,drhow,P,rho,drho,Cx_choice,doNcentered,LZx,Ex)

! Compute the exchange individual energy

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
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Pw(nBas,nBas,nspin)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)
  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin,nEns)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Output variables

  double precision,intent(out)  :: LZx(nspin)
  double precision,intent(out)  :: Ex(nspin,nEns)

  select case (rung)

!   Hartree calculation

    case(0) 

      Ex = 0d0

!   LDA functionals

    case(1) 

      call unrestricted_lda_exchange_individual_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,&
                                                       rhow,rho,Cx_choice,doNcentered,LZx,Ex)

!   GGA functionals

    case(2) 

      call unrestricted_gga_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,LZx,Ex)

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,LZx,Ex)

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,nBas,ERI,Pw,rhow,drhow,P,rho,drho,LZx,Ex)

  end select
 
end subroutine unrestricted_exchange_individual_energy
