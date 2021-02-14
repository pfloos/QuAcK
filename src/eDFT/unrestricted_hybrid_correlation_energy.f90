subroutine unrestricted_hybrid_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute the unrestricted version of the correlation energy for hybrid functionals

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
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

    case('HF') 

      Ec(:) = 0d0

    case('LYP') 

      aC = 0.81d0

      call unrestricted_lda_correlation_energy('VWN5        ',nEns,wEns,nGrid,weight,rho,EcLDA)
      call unrestricted_gga_correlation_energy('LYP         ',nEns,wEns,nGrid,weight,rho,drho,EcGGA)

      Ec(:) = EcLDA(:) + aC*(EcGGA(:) - EcLDA(:)) 

    case('PBE') 

      call unrestricted_gga_correlation_energy('PBE         ',nEns,wEns,nGrid,weight,rho,drho,EcGGA)

      Ec(:) = EcGGA(:)

    case default
  
      call print_warning('!!! Hybrid correlation energy not available !!!')
      stop

  end select
 
end subroutine unrestricted_hybrid_correlation_energy
