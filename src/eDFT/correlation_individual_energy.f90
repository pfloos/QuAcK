subroutine correlation_individual_energy(rung,DFA,nEns,wEns,nGrid,weight,rhow,drhow,rho,drho,Ec)

! Compute the correlation energy of individual states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  double precision              :: EcLDA(nsp)
  double precision              :: EcGGA(nsp)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec(nsp)

  select case (rung)

!   Hartree calculation

    case(0) 

      Ec(:) = 0d0

!   LDA functionals

    case(1) 

      call lda_correlation_individual_energy(DFA,nEns,wEns(:),nGrid,weight(:),rhow(:,:),rho(:,:),Ec(:))

!   GGA functionals

    case(2) 

      call print_warning('!!! Individual energies NYI for GGAs !!!')
      stop

!   Hybrid functionals

    case(4) 

      call print_warning('!!! Individual energies NYI for hybrids !!!')
      stop

      aC = 0.81d0

!   Hartree-Fock calculation

    case(666) 

      Ec(:) = 0d0

  end select
 
end subroutine correlation_individual_energy
