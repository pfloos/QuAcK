function exchange_energy(rung,nGrid,weight,nBas,P,FxHF,rho,drho) result(Ex)

! Compute the exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: FxHF(nBas,nBas)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  double precision              :: ExLDA,ExGGA,ExHF
  double precision              :: cX,aX,aC
  double precision              :: Ex

! Output variables

! Memory allocation

  Ex    = 0d0
  ExLDA = 0d0
  ExGGA = 0d0
  ExHF  = 0d0

  select case (rung)

!   Hartree calculation
    case(0) 

      Ex = 0d0

!   LDA functionals
    case(1) 

      call lda_exchange_energy(nGrid,weight,rho,ExLDA)

      Ex = ExLDA

!   GGA functionals
    case(2) 

      call gga_exchange_energy(nGrid,weight,rho,drho,ExGGA)

      Ex = ExGGA

!   Hybrid functionals
    case(4) 

      cX = 0.20d0
      aX = 0.72d0
      aC = 0.81d0

      call lda_exchange_energy(nGrid,weight,rho,ExLDA)
      call gga_exchange_energy(nGrid,weight,rho,drho,ExGGA)
      call fock_exchange_energy(nBas,P,FxHF,ExHF)

      Ex = ExLDA              & 
         + cX*(ExHF  - ExLDA) &
         + aX*(ExGGA - ExLDA) 

!   Hartree-Fock calculation
    case(666) 

      call fock_exchange_energy(nBas,P,FxHF,ExHF)

      Ex = ExHF

  end select
 
end function exchange_energy
