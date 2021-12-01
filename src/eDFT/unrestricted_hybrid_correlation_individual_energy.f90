subroutine unrestricted_hybrid_correlation_individual_energy(DFA,nEns,wEns,nGrid,weight, & 
                                                             rhow,drhow,rho,drho,doNcentered,LZc,Ec)

! Compute the hybrid correlation energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  double precision,intent(in)   :: rho(nGrid,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nEns)
  logical,intent(in)            :: doNcentered

! Output variables

  double precision              :: LZc(nsp)
  double precision              :: Ec(nsp,nEns)

! Select correlation functional

  select case (DFA)

    case (1)

      LZc(:)  = 0d0
      Ec(:,:) = 0d0

    case default

      call print_warning('!!! Hybrid correlation individual energy not available !!!')
      stop

  end select

end subroutine unrestricted_hybrid_correlation_individual_energy
