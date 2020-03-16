subroutine exchange_Levy_Zahariev_shift(rung,DFA,nEns,wEns,nGrid,weight,rho,drho,ExLZ)

! Compute the exchange part of the Levy-Zahariev shift

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


! Output variables

  double precision,intent(out)  :: ExLZ

  select case (rung)

!   Hartree calculation

    case(0) 

      ExLZ = 0d0

!   LDA functionals

    case(1) 

      call lda_exchange_Levy_Zahariev_shift(DFA,nEns,wEns(:),nGrid,weight(:),rho(:),ExLZ)

!   GGA functionals

    case(2) 

      call print_warning('!!! Exchange LZ shift NYI for GGAs !!!')
      stop

!   Hybrid functionals

    case(4) 

      call print_warning('!!! Exchange LZ shift NYI for hybrids !!!')
      stop

!   Hartree-Fock calculation

    case(666) 

      ExLZ = 0d0

  end select
 
end subroutine exchange_Levy_Zahariev_shift
