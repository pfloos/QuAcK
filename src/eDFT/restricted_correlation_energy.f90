subroutine restricted_correlation_energy(rung,DFA,LDA_centered,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute the correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  double precision              :: EcLDA
  double precision              :: EcGGA
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec

  select case (rung)

!   Hartree calculation

    case(0) 

      Ec = 0d0

!   LDA functionals

    case(1) 

      call restricted_lda_correlation_energy(DFA,LDA_centered,nEns,wEns(:),nGrid,weight(:),rho(:),Ec)

!   GGA functionals

    case(2) 

      call print_warning('!!! restricted correlation energies NYI for GGAs !!!')
      stop


    case(4) 

      call print_warning('!!! restricted correlation energies NYI for Hybrids !!!')
      stop

!   Hartree-Fock calculation

    case(666) 

      Ec = 0d0

  end select
 
end subroutine restricted_correlation_energy
