subroutine hybrid_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute the unrestricted version of the correlation energy for hybrid functionals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  double precision              :: EcLDA(nsp)
  double precision              :: EcGGA(nsp)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec(nsp)

  select case (DFA)

    case(1) 

      Ec(:) = 0d0

    case(2) 

      aC = 0.81d0

      call lda_correlation_energy(3,nEns,wEns,nGrid,weight,rho,EcLDA)
      call gga_correlation_energy(1,nEns,wEns,nGrid,weight,rho,drho,EcGGA)

      Ec(:) = EcLDA(:) + aC*(EcGGA(:) - EcLDA(:)) 

    case(3) 

      call gga_correlation_energy(1,nEns,wEns,nGrid,weight,rho,drho,Ec)

    case(4) 

      call gga_correlation_energy(2,nEns,wEns,nGrid,weight,rho,drho,Ec)

    case default
  
      call print_warning('!!! Hybrid correlation energy not available !!!')
      stop

  end select
 
end subroutine hybrid_correlation_energy
