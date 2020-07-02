subroutine exchange_energy(rung,DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas,P,FxHF,rho,drho,Ex)

! Compute the exchange energy

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
  double precision,intent(in)   :: FxHF(nBas,nBas)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision              :: ExLDA,ExGGA,ExHF
  double precision              :: cX,aX,aC

! Output variables

  double precision,intent(out)  :: Ex

  select case (rung)

!   Hartree calculation

    case(0) 

      Ex = 0d0

!   LDA functionals

    case(1) 

      call lda_exchange_energy(DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rho,ExLDA)

      Ex = ExLDA

!   GGA functionals

    case(2) 

      call gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,ExGGA)

      Ex = ExGGA

!   Hybrid functionals

    case(4) 

      cX = 0.20d0
      aX = 0.72d0
      aC = 0.81d0

      call lda_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,ExLDA)
      call gga_exchange_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,ExGGA)
      call fock_exchange_energy(nBas,P,FxHF,ExHF)

      Ex = ExLDA              & 
         + cX*(ExHF  - ExLDA) &
         + aX*(ExGGA - ExLDA) 

!   Hartree-Fock calculation

    case(666) 

      call fock_exchange_energy(nBas,P,FxHF,ExHF)

      Ex = ExHF

  end select
 
end subroutine exchange_energy
