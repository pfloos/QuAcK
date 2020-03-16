subroutine exchange_individual_energy(rung,DFA,nEns,wEns,nGrid,weight,nBas, & 
                                      ERI,P,FxHF,rhow,drhow,rho,drho,Ex)

! Compute the exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: FxHF(nBas,nBas)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision              :: ExLDA
  double precision              :: ExHF

! Output variables

  double precision,intent(out)  :: Ex

  select case (rung)

!   Hartree calculation

    case(0) 

      Ex = 0d0

!   LDA functionals

    case(1) 

      call lda_exchange_individual_energy(DFA,nEns,wEns,nGrid,weight,rhow,rho,ExLDA)

      Ex = ExLDA

!   GGA functionals

    case(2) 

      call print_warning('!!! Individual energies NYI for GGAs !!!')
      stop

!   Hybrid functionals

    case(4) 

      call print_warning('!!! Individual energies NYI for Hybrids !!!')
      stop

!   Hartree-Fock calculation

    case(666) 

      call fock_exchange_potential(nBas,P,ERI,FxHF)
      call fock_exchange_energy(nBas,P,FxHF,ExHF)

      Ex = ExHF

  end select
 
end subroutine exchange_individual_energy
