subroutine restricted_correlation_Levy_Zahariev_shift(rung,DFA,nEns,wEns,nGrid,weight,rho,drho,EcLZ)

! Compute the correlation part of the Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision              :: EcLZLDA,EcLZGGA
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: EcLZ

  select case (rung)

!   Hartree calculation

    case(0) 

      EcLZ = 0d0

!   LDA functionals

    case(1) 

      call restricted_lda_correlation_Levy_Zahariev_shift(DFA,nEns,wEns(:),nGrid,weight(:),rho(:),EcLZ)

!   GGA functionals

    case(2) 

      call print_warning('!!! Individual energies NYI for GGAs !!!')
      stop

!   Hybrid functionals

    case(4) 

      call print_warning('!!! Individual energies NYI for hybrids !!!')
      stop

      aC = 0.81d0

      EcLZ = EcLZLDA + aC*(EcLZGGA - EcLZLDA) 

!   Hartree-Fock calculation

    case(666) 

      EcLZ = 0d0

  end select
 
end subroutine restricted_correlation_Levy_Zahariev_shift
